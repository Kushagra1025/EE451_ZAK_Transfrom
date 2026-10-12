%% S8_full_chain.m -- speech -> PCM -> 3GPP 16QAM -> DD -> IZak -> AWGN -> Zak -> EQ -> demap -> speech
clear; close all; clc;
rng(7);

M = 32; N = 16; D = 4; Lg = 65;  K = M*N;          % grid, oversampling, taps, symbols/frame
fprintf('M=%d N=%d D=%d Lg=%d  symbols/frame=%d  bits/frame=%d\n', M, N, D, Lg, K, 4*K);
P = make_pulses(M, N, D, Lg);
Qf = @(x) 0.5*erfc(x/sqrt(2));
ber_th = @(ebn0_db) arrayfun(@(e) ...
    0.75*Qf(sqrt(0.8*10^(e/10))) + 0.5*Qf(3*sqrt(0.8*10^(e/10))) - 0.25*Qf(5*sqrt(0.8*10^(e/10))), ebn0_db);

%% ---------- Part A: 3GPP 16QAM constellation ----------
allb = zeros(4,16);
for i = 0:15, allb(:,i+1) = bitget(i, 4:-1:1).'; end        % b0 b1 b2 b3 (b0 = MSB of the label)
cpt  = qam16_map(allb);
fprintf('\n--- Part A: 16QAM (TS 38.211 clause 5.1.4 as implemented here) ---\n');
fprintf('distinct points         : %d\n', numel(unique(round(cpt*sqrt(10)*1e6))));
fprintf('mean energy             : %.4f\n', mean(abs(cpt).^2));
fprintf('demap(map(b)) == b      : %d\n', isequal(qam16_demap(cpt), allb));
hd = [];
for i = 1:16, for j = 1:16
    if i ~= j && abs(abs(cpt(i)-cpt(j))*sqrt(10) - 2) < 1e-9, hd(end+1) = sum(allb(:,i) ~= allb(:,j)); end %#ok<SAGROW>
end, end
fprintf('Hamming distance between nearest neighbours: %s\n', mat2str(unique(hd)));
fprintf('label b0b1b2b3 -> point*sqrt(10)\n');
for i = 1:16
    fprintf('   %d%d%d%d -> %+d%+dj\n', allb(:,i), round(real(cpt(i))*sqrt(10)), round(imag(cpt(i))*sqrt(10)));
end
b6 = [1 0 1 0 0 1];            % slide p14 example, 64QAM
I6 = (1-2*b6(1))*(4-(1-2*b6(3))*(2-(1-2*b6(5))));
Q6 = (1-2*b6(2))*(4-(1-2*b6(4))*(2-(1-2*b6(6))));
fprintf('64QAM bits 101001 by the 5.1.5 formula (before 1/sqrt(42)): %+d%+dj   (slide p14 shows 5+3j)\n', I6, Q6);
if exist('nrSymbolModulate','file') ~= 0
    try
        bb = randi([0 1], 4*1000, 1);
        d1 = nrSymbolModulate(bb, '16QAM');
        d2 = qam16_map(reshape(bb,4,[]));
        fprintf('cross-check with nrSymbolModulate: max diff = %.2e\n', max(abs(d1(:) - d2(:))));
    catch
        fprintf('nrSymbolModulate present but not usable (license?)\n');
    end
end
figure('Name','Fig1 16QAM constellation'); hold on; grid on;
plot(real(cpt), imag(cpt), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
for i = 1:16
    text(real(cpt(i))+0.03, imag(cpt(i))+0.05, sprintf('%d%d%d%d', allb(:,i)), 'FontSize', 9);
end
axis equal; xlim([-1.3 1.3]); ylim([-1.3 1.3]);
xlabel('I'); ylabel('Q'); title('3GPP 16QAM labels b_0b_1b_2b_3');

%% ---------- Part B: speech -> PCM -> chain ----------
fs = 8000;  Ns = 8192;  n = (0:Ns-1).';  t = n/fs;
f0 = 120 + 60*t/1.024;  ph = 2*pi*cumsum(f0)/fs;
x = zeros(Ns,1);
for k = 1:21, x = x + sin(k*ph)/k; end
x = x .* (0.55 + 0.45*sin(2*pi*3*t - pi/2));
x = 0.9*x/max(abs(x));
qi = min(max(round((x+1)/2*255), 0), 255);                 % 8-bit PCM
bitsM = zeros(8, Ns);
for k = 1:8, bitsM(k,:) = bitget(qi.', 9-k); end            % MSB first
bits = bitsM(:);                                            % 65536 bits
B4 = reshape(bits, 4, []);                                  % 4 x (K*32)
nfr = size(B4,2)/K;
xr0 = qi/255*2 - 1;
sqnr = @(a,b) 10*log10(mean((a-mean(a)).^2)/mean((a-b).^2));
fprintf('\n--- Part B: speech/PCM through the chain ---\n');
fprintf('samples=%d  bits=%d  frames=%d  signal rms=%.4f\n', Ns, numel(bits), nfr, std(x,1));
fprintf('8-bit PCM SQNR (error-free)       : %.2f dB\n', sqnr(x, xr0));

% Level 0: no noise, rect pulse
bh = zeros(size(B4));
for f = 1:nfr
    cols = (f-1)*K + (1:K);
    bh(:,cols) = chain_frame(B4(:,cols), P.rect, M, N, D, 0, []);
end
fprintf('Level-0 loopback (noise-free) BER : %g   bit-exact: %d\n', mean(bh(:) ~= bits), isequal(bh(:), bits));

ebv = [4 8 12 16];  nrep = 10;
fprintf('\n Eb/N0   BER(sim)   BER(theory)   SQNR mean   min    max   [dB]\n');
xr_keep = zeros(Ns, numel(ebv));
for ie = 1:numel(ebv)
    es  = 4*10^(ebv(ie)/10);  st2 = 1/(N*es);
    bers = zeros(nrep,1);  sq = zeros(nrep,1);
    for r = 1:nrep
        bh = zeros(size(B4));
        for f = 1:nfr
            cols = (f-1)*K + (1:K);
            bh(:,cols) = chain_frame(B4(:,cols), P.rect, M, N, D, st2, []);
        end
        bers(r) = mean(bh(:) ~= bits);
        qh = (2.^(7:-1:0)) * reshape(bh(:), 8, []);
        xr = qh.'/255*2 - 1;
        sq(r) = sqnr(x, xr);
        if r == 1, xr_keep(:,ie) = xr; end
    end
    fprintf('%6.0f  %9.2e  %11.2e   %8.2f  %6.2f %6.2f\n', ebv(ie), mean(bers), ber_th(ebv(ie)), mean(sq), min(sq), max(sq));
end

% Fig 2: PCM quantisation
figure('Name','Fig2 speech and PCM');
subplot(3,1,1); plot(t, x); grid on; ylabel('x'); title('synthetic speech-like signal, 8 kHz');
i0 = 1200:1440;
subplot(3,1,2); plot(t(i0), x(i0), '-', t(i0), xr0(i0), '.-'); grid on; legend({'analog','8-bit PCM'}); ylabel('zoom (30 ms)');
subplot(3,1,3); plot(t, x - xr0); grid on; ylabel('quantisation error'); xlabel('time (s)');

% Fig 3: rect constellations after the Zak (8 frames each)
figure('Name','Fig3 received DD symbols (rect)');
for ie = 1:numel(ebv)
    es = 4*10^(ebv(ie)/10);  st2 = 1/(N*es);  yy = [];
    for f = 1:8
        cols = (f-1)*K + (1:K);
        [~, y] = chain_frame(B4(:,cols), P.rect, M, N, D, st2, []);
        yy = [yy; y]; %#ok<AGROW>
    end
    subplot(2,2,ie); plot(real(yy), imag(yy), '.', 'MarkerSize', 4); hold on; grid on;
    plot(real(cpt), imag(cpt), 'wo', 'MarkerSize', 6, 'LineWidth', 1);
    axis equal; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
    title(sprintf('Eb/N0 = %d dB', ebv(ie))); xlabel('I'); ylabel('Q');
end

% Fig 4: reconstructed speech
figure('Name','Fig4 reconstructed speech');
iz = 1:2000;
for ie = 1:numel(ebv)
    subplot(4,1,ie); plot(t(iz), x(iz), '-', t(iz), xr_keep(iz,ie), '-'); grid on;
    ylabel(sprintf('%d dB', ebv(ie))); ylim([-1.2 1.2]);
    if ie == 1, title('original vs received (first 0.25 s)'); legend({'original','received'}); end
end
xlabel('time (s)');

%% ---------- Part C: DD equaliser and BER sweep ----------
fprintf('\n--- Part C: H matrices (S4 theory) ---\n');
H_gs = build_H(P.gs.g,   P.gs.ell0,   D, M, N);
H_or = build_H(P.orth.g, P.orth.ell0, D, M, N);
Xt = reshape(qam16_map(randi([0 1],4,K)), M, N);
Yt = zak_sampled(izak_sampled(Xt, P.gs.g, P.gs.ell0, D), M, N, P.gs.g, P.gs.ell0, D);
fprintf('H*X vs chain output (gs, noise-free): %.2e\n', max(abs(Yt(:) - H_gs*Xt(:))));
for nm = {'gs','orth'}
    Hh = H_gs; if strcmp(nm{1},'orth'), Hh = H_or; end
    ev = eig((Hh+Hh')/2);
    fprintf('%-5s  eig(H): min %.4f max %.4f  mean(1/eig) %.4f  ||H-I||_F/sqrt(K) = %.4f\n', ...
        nm{1}, min(ev), max(ev), mean(1./ev), norm(Hh - eye(K),'fro')/sqrt(K));
end

ebs = 0:2:12;
names = {'rect','gs','gs+MMSE','orth'};
ber = zeros(numel(ebs), 4);  nbt = zeros(numel(ebs), 4);
for ie = 1:numel(ebs)
    es = 4*10^(ebs(ie)/10);  st2 = 1/(N*es);
    Wm = inv(H_gs + (1/es)*eye(K));                       % MMSE for noise covariance (1/es)*H
    for c = 1:4
        switch c
            case 1, Pc = P.rect; W = [];
            case 2, Pc = P.gs;   W = [];
            case 3, Pc = P.gs;   W = Wm;
            case 4, Pc = P.orth; W = [];
        end
        nerr = 0; nb = 0;
        for f = 1:600
            b = randi([0 1], 4, K);
            bh = chain_frame(b, Pc, M, N, D, st2, W);
            nerr = nerr + sum(bh(:) ~= b(:));  nb = nb + numel(b);
            if nerr >= 200 && f >= 20, break; end
        end
        ber(ie,c) = nerr/nb;  nbt(ie,c) = nb;
    end
end
fprintf('\nBER vs Eb/N0 (stops at 200 errors or 600 frames)\n');
fprintf('Eb/N0    theory      rect       gs      gs+MMSE     orth\n');
for ie = 1:numel(ebs)
    fprintf('%5.0f  %9.2e  %9.2e  %9.2e  %9.2e  %9.2e\n', ebs(ie), ber_th(ebs(ie)), ber(ie,:));
end
figure('Name','Fig5 BER vs Eb/N0'); hold on; grid on;
eg = 0:0.25:14;
semilogy(eg, ber_th(eg), 'w-', 'LineWidth', 1.5);
mk = {'o-','s-','^-','d-'};
for c = 1:4
    yv = ber(:,c);  yv(yv == 0) = NaN;
    semilogy(ebs, yv, mk{c}, 'LineWidth', 1.4, 'MarkerSize', 7);
end
set(gca,'YScale','log'); ylim([1e-5 1]); xlim([0 14]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend([{'16QAM Gray theory'}, names], 'Location','southwest');
title('End-to-end BER, OTFS chain over AWGN');

% Fig 6: equalised constellations at 12 dB
es = 4*10^(12/10);  st2 = 1/(N*es);  Wm = inv(H_gs + (1/es)*eye(K));
cfg = {P.gs, P.gs, P.orth};  Wc = {[], Wm, []};  ttl = {'gs, no EQ','gs + MMSE','orth, no EQ'};
figure('Name','Fig6 equalised constellations 12 dB');
for c = 1:3
    yy = [];
    for f = 1:8
        [~, y] = chain_frame(randi([0 1],4,K), cfg{c}, M, N, D, st2, Wc{c});
        yy = [yy; y]; %#ok<AGROW>
    end
    subplot(1,3,c); plot(real(yy), imag(yy), '.', 'MarkerSize', 4); hold on; grid on;
    plot(real(cpt), imag(cpt), 'wo', 'MarkerSize', 6, 'LineWidth', 1);
    axis equal; xlim([-1.6 1.6]); ylim([-1.6 1.6]); title(ttl{c}); xlabel('I'); ylabel('Q');
end

%% ---------- Part D: noise-free floor and eigenvalues ----------
fprintf('\n--- Part D: noise-free BER (500 frames = 1.02 Mbit) ---\n');
Wzf = inv(H_gs);
cases = {'rect',P.rect,[]; 'gs, no EQ',P.gs,[]; 'gs, ZF',P.gs,Wzf; 'orth, no EQ',P.orth,[]};
for c = 1:4
    nerr = 0;
    for f = 1:500
        b = randi([0 1],4,K);
        bh = chain_frame(b, cases{c,2}, M, N, D, 0, cases{c,3});
        nerr = nerr + sum(bh(:) ~= b(:));
    end
    fprintf('%-12s  errors = %3d   BER = %.2e\n', cases{c,1}, nerr, nerr/(500*4*K));
end
figure('Name','Fig7 eigenvalues of H'); hold on; grid on;
plot(sort(eig((H_gs+H_gs')/2)), '.', 'MarkerSize', 8);
plot(sort(eig((H_or+H_or')/2)), '.', 'MarkerSize', 8);
yline(1, '--', 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
xlabel('index'); ylabel('eigenvalue of H'); legend({'Gauss*sinc','orthogonalised'}, 'Location','southeast');
title('Eigenvalues of H = Zak o IZak  (all = 1 when the pulse is orthogonal)');

%% ======================= local functions =======================
function d = qam16_map(b)                     % b: 4 x K, rows b0 b1 b2 b3
    d = ((1-2*b(1,:)).*(2-(1-2*b(3,:))) + 1j*(1-2*b(2,:)).*(2-(1-2*b(4,:)))) / sqrt(10);
    d = d(:);
end

function b = qam16_demap(y)
    y = y(:);
    re = real(y)*sqrt(10);  im = imag(y)*sqrt(10);
    b = double([(re < 0).'; (im < 0).'; (abs(re) > 2).'; (abs(im) > 2).']);
end

function [bh, y] = chain_frame(b, P1, M, N, D, st2, W)
    X = reshape(qam16_map(b), M, N);
    s = izak_sampled(X, P1.g, P1.ell0, D);
    w = sqrt(st2/2) * (randn(size(s)) + 1j*randn(size(s)));
    Y = zak_sampled(s + w, M, N, P1.g, P1.ell0, D);
    y = Y(:);
    if ~isempty(W), y = W*y; end
    bh = qam16_demap(y);
end

function H = build_H(g, ell0, D, M, N)
    K = M*N;  H = zeros(K, K);
    for i = 1:K
        X = zeros(M, N);  X(i) = 1;
        Y = zak_sampled(izak_sampled(X, g, ell0, D), M, N, g, ell0, D);
        H(:,i) = Y(:);
    end
end

function P = make_pulses(M, N, D, Lg)
    L = M*N*D;  c = (Lg-1)/2;
    P.rect.g = ones(D,1)/sqrt(D);  P.rect.ell0 = 0;
    t = (-c:c).'/D;
    g = exp(-pi*4*(t/M).^2) .* sincn(t);
    P.gs.g = g/norm(g);  P.gs.ell0 = -c;
    Lb = 2*floor(min(32*D, L-1)/2) + 1;  cb = (Lb-1)/2;
    tb = (-cb:cb).'/D;
    base = exp(-pi*16*(tb/M).^2) .* sincn(tb);  base = base/norm(base);
    go = orth_full(base, D, L);
    gg = go(mod(-c:c, L) + 1);
    P.orth.g = gg(:)/norm(gg);  P.orth.ell0 = -c;
end

function go = orth_full(base, D, L)
    n = numel(base);  c = (n-1)/2;
    full = zeros(L,1);  full(mod(-c:c, L) + 1) = base;
    G = fft(full);  Lq = L/D;  Fo = zeros(Lq,1);
    for k = 0:D-1, Fo = Fo + abs(G(k*Lq + (1:Lq))).^2; end
    go = real(ifft(G .* sqrt(D ./ repmat(Fo, D, 1))));
end

function y = sincn(t)
    y = ones(size(t));  nz = (t ~= 0);
    y(nz) = sin(pi*t(nz)) ./ (pi*t(nz));
end

function s = izak_sampled(X, g, ell0, D)
    [M, N] = size(X);  P = M*D;  L = N*P;
    V = ifft(X, [], 2);
    up = zeros(L, 1);
    k = (0:M-1).' + M*(0:N-1);
    up(D*k(:) + 1) = V(:);
    s = zeros(L, 1);
    for i = 1:numel(g)
        s = s + g(i) * circshift(up, ell0 + i - 1);
    end
end

function Y = zak_sampled(r, M, N, g, ell0, D)
    L = M*N*D;  r = r(:);
    acc = zeros(L, 1);
    for i = 1:numel(g)
        acc = acc + circshift(r, -(ell0 + i - 1)) * conj(g(i));
    end
    k = (0:M-1).' + M*(0:N-1);
    U = reshape(acc(D*k(:) + 1), M, N);
    Y = fft(U, [], 2);
end