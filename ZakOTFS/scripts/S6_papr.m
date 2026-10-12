%% S6_papr.m  -- PAPR of OTFS (rect / Gauss*sinc / orthogonalised) vs OFDM
clear; close all; clc;
rng(1);

M = 32;  N = 16;  D = 4;  Lg = 65;          % grid, oversampling, pulse taps
F = 1000;                                   % random frames for the statistics
L = M*N*D;                                  % samples per frame
fprintf('M=%d N=%d D=%d  samples/frame=%d  frames=%d\n', M, N, D, L, F);

%% ---------- Part A: PAPR statistics over random 16QAM frames ----------
[P.rect, P.gs, P.orth] = deal([]);
P = make_pulses(M, N, D, Lg);
names = {'rect','gs','orth'};
papr_v = struct('rect',zeros(F,1),'gs',zeros(F,1),'orth',zeros(F,1), ...
                'ofdm',zeros(F,1),'ofdm_sym',zeros(F,1));
for f = 1:F
    X = qam16([M N]);
    for k = 1:3
        nm = names{k};
        s = izak_sampled(X, P.(nm).g, P.(nm).ell0, D);
        papr_v.(nm)(f) = papr_db(s);
    end
    so = ofdm_tx(X, D);
    papr_v.ofdm(f) = papr_db(so);
    ps = zeros(N,1);
    for i = 1:N, ps(i) = papr_db(so((i-1)*M*D + (1:M*D))); end
    papr_v.ofdm_sym(f) = mean(ps);
end
fn = {'rect','gs','orth','ofdm','ofdm_sym'};
fprintf('\n%-9s %8s %8s %8s %8s   (dB)\n','waveform','mean','99%','99.9%','max');
for k = 1:numel(fn)
    v = papr_v.(fn{k});
    fprintf('%-9s %8.2f %8.2f %8.2f %8.2f\n', fn{k}, mean(v), prctile_(v,99), prctile_(v,99.9), max(v));
end

%% Fig 1: one frame, instantaneous power
rng(5);  X1 = qam16([M N]);
s_r = izak_sampled(X1, P.rect.g, P.rect.ell0, D);
s_o = izak_sampled(X1, P.orth.g, P.orth.ell0, D);
s_f = ofdm_tx(X1, D);
sigs = {s_r, s_o, s_f};  ttl = {'OTFS rect','OTFS orthogonalised','OFDM'};
figure('Name','Fig1 waveform power');
for k = 1:3
    subplot(3,1,k);
    p = abs(sigs{k}).^2;  tt = (0:L-1)/D;
    plot(tt, 10*log10(p/mean(p)), 'LineWidth', 0.6); hold on; grid on;
    yline(0,'k--'); yline(papr_db(sigs{k}),'r--');
    title(sprintf('%s: PAPR = %.2f dB', ttl{k}, papr_db(sigs{k})));
    ylabel('power / mean (dB)'); ylim([-30 14]);
    if k==3, xlabel('time (Ts)'); end
end

%% Fig 2: CCDF
figure('Name','Fig2 CCDF'); hold on; grid on;
cols = lines(5);
for k = 1:numel(fn)
    x = sort(papr_v.(fn{k}));
    y = (F:-1:1)'/F;
    stairs(x, y, 'LineWidth', 1.6, 'Color', cols(k,:));
end
set(gca,'YScale','log'); ylim([1/F 1]);
xlabel('PAPR_0 (dB)'); ylabel('Pr(PAPR > PAPR_0)');
legend({'OTFS rect','OTFS Gauss*sinc','OTFS orthogonalised','OFDM (whole frame)','OFDM (per symbol, mean)'}, 'Location','southwest');
title('CCDF of PAPR per frame');

%% Fig 3: histograms
figure('Name','Fig3 PAPR histograms'); hold on; grid on;
edges = 4:0.25:12;
for k = 1:4
    histogram(papr_v.(fn{k}), edges, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.5, 'EdgeColor',cols(k,:));
end
xlabel('PAPR (dB)'); ylabel('pdf');
legend({'rect','Gauss*sinc','orthogonalised','OFDM'});
title('PAPR distribution over frames');

%% ---------- Part B: how many DD points are occupied ----------
fprintf('\nSingle DD point: expected 10*log10(M) = %.2f dB\n', 10*log10(M));
X = zeros(M,N);  X(6,4) = 1;
for k = 1:3
    nm = names{k};
    s = izak_sampled(X, P.(nm).g, P.(nm).ell0, D);
    fprintf('  %-5s single point PAPR = %.2f dB\n', nm, papr_db(s));
end
% Fig 4: waveform of single point (rect)
s1 = izak_sampled(X, P.rect.g, P.rect.ell0, D);
figure('Name','Fig4 single DD point');
stem((0:3*M*D-1)/D, abs(s1(1:3*M*D)).^2, 'Marker','none'); grid on;
xlabel('time (Ts)'); ylabel('|s|^2');
title('One DD point: one pulse per Zak period (spacing M*Ts), duty cycle 1/M');

Ks = [1 2 4 8 32 128 512];  nrep = 100;
fill_p = zeros(numel(Ks), 2);
for ik = 1:numel(Ks)
    K = Ks(ik);
    for kk = 1:2
        nm = names{1 + 2*(kk-1)};            % 'rect' then 'orth'
        acc = 0;
        for r = 1:nrep
            Xv = zeros(M*N,1);  idx = randperm(M*N, K);
            Xv(idx) = qam16([K 1]);
            s = izak_sampled(reshape(Xv,M,N), P.(nm).g, P.(nm).ell0, D);
            acc = acc + papr_db(s);
        end
        fill_p(ik,kk) = acc/nrep;
    end
end
fprintf('\nFill sweep (mean PAPR dB over %d frames)\n   K     rect    orth\n', nrep);
fprintf('%4d  %6.2f  %6.2f\n', [Ks(:) fill_p]');
figure('Name','Fig5 PAPR vs occupied DD points');
semilogx(Ks, fill_p(:,1), 'o-', Ks, fill_p(:,2), 's-', 'LineWidth', 1.6); grid on; hold on;
yline(10*log10(M), '--', '10log_{10}M', 'Color',[0.7 0.7 0.7], 'HandleVisibility','off');
xlabel('number of occupied DD points K (out of MN)'); ylabel('mean PAPR (dB)');
legend({'rect','orthogonalised'},'Location','northeast');
title('Sparse DD occupancy = high PAPR');

%% ---------- Part C: sweep N (Doppler bins) ----------
Ns = [1 2 4 8 16 32];  nrep = 200;
Pn = zeros(numel(Ns), 4);
for in = 1:numel(Ns)
    Nn = Ns(in);  PP = make_pulses(M, Nn, D, Lg);
    for r = 1:nrep
        Xn = qam16([M Nn]);
        for k = 1:3
            nm = names{k};
            Pn(in,k) = Pn(in,k) + papr_db(izak_sampled(Xn, PP.(nm).g, PP.(nm).ell0, D));
        end
        Pn(in,4) = Pn(in,4) + papr_db(ofdm_tx(Xn, D));
    end
end
Pn = Pn/nrep;
fprintf('\nN sweep (M=%d, D=%d, mean PAPR dB, %d frames)\n    N    rect      gs    orth    ofdm\n', M, D, nrep);
fprintf('%4d  %6.2f  %6.2f  %6.2f  %6.2f\n', [Ns(:) Pn]');
figure('Name','Fig6 PAPR vs N');
semilogx(Ns, Pn, 'o-', 'LineWidth', 1.6); grid on;
xlabel('N (Doppler bins / Zak periods)'); ylabel('mean PAPR (dB)');
legend({'rect','Gauss*sinc','orthogonalised','OFDM'},'Location','southeast');
title('N = 1 is a single-carrier stream; PAPR grows with N');

%% ---------- Part D: spectrum (what the pulse actually buys) ----------
fprintf('\nSpectrum (0 dB = mean level for |f|<0.4/Ts)\n');
fax = ((0:L-1) - L/2) / L * D;               % frequency in units of 1/Ts
figure('Name','Fig7 PSD'); hold on; grid on;
for k = 1:3
    nm = names{k};  S = zeros(L,1);
    for r = 1:200
        s = izak_sampled(qam16([M N]), P.(nm).g, P.(nm).ell0, D);
        S = S + abs(fft(s)).^2;
    end
    S = fftshift(S/200);
    S = S / mean(S(abs(fax) < 0.4));
    plot(fax, 10*log10(S + eps), 'LineWidth', 1.3);
    at = @(x) 10*log10(S(find(abs(fax-x) == min(abs(fax-x)), 1)));
    fprintf('  %-5s  f=0.75/Ts: %6.1f dB   f=1.5/Ts: %6.1f dB\n', nm, at(0.75), at(1.5));
end
xlabel('frequency (1/Ts)'); ylabel('PSD (dB)'); ylim([-70 5]); xlim([-D/2 D/2]);
xline(-0.5,'k:'); xline(0.5,'k:');
legend({'rect','Gauss*sinc','orthogonalised'});
title('Averaged PSD; band edge at +/-0.5/Ts');

%% ======================= local functions =======================
function X = qam16(sz)
    X = ((2*randi([0 3],sz)-3) + 1j*(2*randi([0 3],sz)-3)) / sqrt(10);
end

function p = papr_db(s)
    a = abs(s(:)).^2;
    p = 10*log10(max(a)/mean(a));
end

function v = prctile_(x, p)          % no toolbox needed
    x = sort(x(:));  n = numel(x);
    v = x(min(n, max(1, ceil(p/100*n))));
end

function s = ofdm_tx(X, D)           % N OFDM symbols, M subcarriers, oversampled by D
    [M, N] = size(X);
    Z = zeros(M*D, N);  Z(1:M,:) = X;
    s = reshape(ifft(Z, [], 1), [], 1);
end

function P = make_pulses(M, N, D, Lg)
    L = M*N*D;  c = (Lg-1)/2;
    P.rect.g = ones(D,1)/sqrt(D);  P.rect.ell0 = 0;
    t = (-c:c).'/D;
    g = exp(-pi*4*(t/M).^2) .* sincn(t);
    P.gs.g = g/norm(g);  P.gs.ell0 = -c;
    Lb = 2*floor(min(32*D, L-1)/2) + 1;  cb = (Lb-1)/2;
    tb = (-cb:cb).'/D;
    base = exp(-pi*16*(tb/M).^2) .* sincn(tb);
    base = base/norm(base);
    go = orth_full(base, D, L);
    gg = go(mod(-c:c, L) + 1);
    P.orth.g = gg(:)/norm(gg);  P.orth.ell0 = -c;
end

function go = orth_full(base, D, L)
    n = numel(base);  c = (n-1)/2;
    full = zeros(L,1);
    full(mod(-c:c, L) + 1) = base;
    G = fft(full);  Lq = L/D;
    Fo = zeros(Lq,1);
    for k = 0:D-1
        Fo = Fo + abs(G(k*Lq + (1:Lq))).^2;
    end
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