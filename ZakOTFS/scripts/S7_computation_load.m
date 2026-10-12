%% S7_computation_load.m  -- operation counts of IZak / Zak (slides p7-8, p16-20)
clear; close all; clc;
rng(2);

%% ---------- Part A: reproduce the slide example ----------
M = 128; N = 128; B = 500e6; Lg = 31;  MN = M*N;
C_dir   = M*N^2;                 % direct IZak MACs / frame
C_fft   = M*(N/2)*log2(N);       % IZak FFT butterflies / frame
C_pulse = MN*Lg;                 % pulse synthesis MACs / frame
C_txfft = C_fft + C_pulse;
C_txdir = C_dir + C_pulse;
f_frame = B/MN;
fprintf('--- Part A: slide example M=N=%d, B=%g MHz, Lg=%d ---\n', M, B/1e6, Lg);
chk('MN',                 MN,        16384);
chk('direct IZak MAC',    C_dir,     2097152);
chk('FFT IZak butterfly', C_fft,     57344);
chk('pulse synthesis',    C_pulse,   507904);
chk('Tx FFT total',       C_txfft,   565248);
chk('frame rate [1/s]',   f_frame,   30517.6);
chk('R IZak FFT [G/s]',   f_frame*C_fft/1e9,   1.75);
chk('R pulse [G MAC/s]',  f_frame*C_pulse/1e9, 15.5);
chk('R Tx direct [G/s]',  B*(N+Lg)/1e9,        79.5);
fprintf('Tx direct / Tx FFT (op count)     = %.2f\n', C_txdir/C_txfft);
fprintf('transform only: direct / FFT      = %.2f\n', C_dir/C_fft);
fprintf('pulse share of FFT-based Tx cost  = %.1f %%\n', 100*C_pulse/C_txfft);
fprintf('Rx: projection mult = %d, FFT twiddle mult = %d, adds = %d\n', ...
        MN*Lg, C_fft, MN*(Lg-1) + MN*log2(N));

%% ---------- Part B: count a real radix-2 FFT ----------
fprintf('\n--- Part B: counted radix-2 DIT FFT ---\n');
fprintf('%5s %9s %9s %12s %14s %10s\n','N','butterf','(N/2)log2N','twiddle!=triv','closed form','max err');
Nlist = 2.^(1:10);
nt_all = zeros(size(Nlist));  nb_all = zeros(size(Nlist));
for i = 1:numel(Nlist)
    Nn = Nlist(i);
    x = randn(Nn,1) + 1j*randn(Nn,1);
    [a, nb, nt] = fft_count(x);
    nt_all(i) = nt;  nb_all(i) = nb;
    cf = (Nn/2)*log2(Nn) - 3*Nn/2 + 2;
    fprintf('%5d %9d %9d %12d %14d %10.1e\n', Nn, nb, (Nn/2)*log2(Nn), nt, cf, max(abs(a - fft(x))));
end
fprintf('N=128: counted nontrivial twiddles per row = %d, per frame (M=128) = %d vs slide bound %d\n', ...
        nt_all(Nlist==128), 128*nt_all(Nlist==128), 57344);

%% ---------- Part C: counted Tx / Rx chain on a small grid ----------
Ms = 16; Ns = 8; Ds = 4; Lgs = 33; Ps = Ms*Ds; Ls = Ns*Ps;
c = (Lgs-1)/2;  t = (-c:c).'/Ds;
g = exp(-pi*4*(t/Ms).^2) .* sincn(t);  g = g/norm(g);  ell0 = -c;
X = qam16([Ms Ns]);
[s, ctx] = tx_counted(X, g, ell0, Ds);
[Y, crx] = rx_counted(s, Ms, Ns, g, ell0, Ds);
s_ref = izak_sampled(X, g, ell0, Ds);
Y_ref = zak_sampled(s, Ms, Ns, g, ell0, Ds);
fprintf('\n--- Part C: counted chain, M=%d N=%d D=%d Lg=%d ---\n', Ms, Ns, Ds, Lgs);
fprintf('Tx output vs izak_sampled: %.2e\n', max(abs(s - s_ref)));
fprintf('Rx output vs zak_sampled : %.2e\n', max(abs(Y(:) - Y_ref(:))));
fprintf('                      counted   formula\n');
fprintf('Tx FFT butterflies  %9d %9d\n', ctx.but, Ms*(Ns/2)*log2(Ns));
fprintf('Tx nontriv twiddles %9d %9s\n', ctx.tw, '-');
fprintf('Tx pulse MACs       %9d %9d\n', ctx.mac, Ms*Ns*Lgs);
fprintf('Rx projection mult  %9d %9d\n', crx.mul, Ms*Ns*Lgs);
fprintf('Rx projection add   %9d %9d\n', crx.add, Ms*Ns*(Lgs-1));
fprintf('Rx FFT butterflies  %9d %9d\n', crx.but, Ms*(Ns/2)*log2(Ns));
fprintf('Rx FFT adds (2/butterfly) %3d %9d\n', 2*crx.but, Ms*Ns*log2(Ns));

%% ---------- Part D: sweeps ----------
% Fig 1: total cost per frame vs N (M=128, Lg=31)
Nv = 2.^(2:10);  Mf = 128;  Lg0 = 31;
Cd = Mf*Nv.^2;  Cf = Mf*(Nv/2).*log2(Nv);  Cp = Mf*Nv*Lg0;
figure('Name','Fig1 cost per frame vs N');
loglog(Nv, Cd+Cp, 'o-', Nv, Cf+Cp, 's-', Nv, Cd, '--', Nv, Cf, '--', Nv, Cp, ':', 'LineWidth', 1.6); grid on;
xlabel('N'); ylabel('operations per frame');
legend({'direct total','FFT total','direct transform only','FFT transform only','pulse synthesis'},'Location','northwest');
title(sprintf('Tx cost per frame, M=%d, L_g=%d', Mf, Lg0));

% Fig 2: per symbol cost and pulse share
figure('Name','Fig2 per-symbol cost and pulse share');
subplot(1,2,1);
loglog(Nv, Nv+Lg0, 'o-', Nv, log2(Nv)/2+Lg0, 's-', Nv, Nv, '--', Nv, log2(Nv)/2, '--', 'LineWidth', 1.6); grid on;
xlabel('N'); ylabel('operations per DD symbol');
legend({'direct (N+L_g)','FFT (log_2N/2+L_g)','direct transform','FFT transform'},'Location','northwest');
title('Cost per symbol');
subplot(1,2,2); hold on; grid on;
for Lgx = [7 15 31 63 127]
    plot(Nv, 100*Lgx./(Lgx + log2(Nv)/2), 'o-', 'LineWidth', 1.4, 'DisplayName', sprintf('L_g=%d', Lgx));
end
set(gca,'XScale','log'); xlabel('N'); ylabel('pulse share of FFT-based Tx cost (%)'); ylim([0 100]);
legend('Location','southeast'); title('Pulse synthesis dominates');

% Fig 3: rate vs D at fixed pulse span (Lg = 2*K*D+1)
Dv = 1:16;  Kv = [4 8 16];  Nfix = 128;
figure('Name','Fig3 rate vs D'); hold on; grid on;
for K = Kv
    plot(Dv, B*(2*K*Dv+1)/1e9, 'o-', 'LineWidth', 1.5, 'DisplayName', sprintf('pulse, span +/-%d Ts', K));
end
yline(B*log2(Nfix)/2/1e9, '--', 'FFT butterflies (independent of D)', 'Color',[0.7 0.7 0.7], 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
plot(4, B*31/1e9, 'p', 'MarkerSize', 14, 'LineWidth', 1.5, 'DisplayName', 'slide: L_g=31');
xlabel('oversampling D'); ylabel('rate (G ops/s), B = 500 MHz'); legend('Location','northwest');
title('D matters only through L_g');

% Fig 4: accuracy eps vs Lg (S4/S6 pulses), M=32, N=16, D=4
Mq = 32; Nq = 16; Dq = 4; Lq = Mq*Nq*Dq;
Lgv = 5:2:257;  e_gs = zeros(size(Lgv));  e_or = zeros(size(Lgv));
Lb = 2*floor(min(32*Dq, Lq-1)/2) + 1;  cb = (Lb-1)/2;  tb = (-cb:cb).'/Dq;
base = exp(-pi*16*(tb/Mq).^2) .* sincn(tb);  base = base/norm(base);
go = orth_full(base, Dq, Lq);
for i = 1:numel(Lgv)
    cc = (Lgv(i)-1)/2;  tt = (-cc:cc).'/Dq;
    g1 = exp(-pi*4*(tt/Mq).^2) .* sincn(tt);  g1 = g1/norm(g1);
    g2 = go(mod(-cc:cc, Lq) + 1);  g2 = g2(:)/norm(g2);
    e_gs(i) = eps_pred(g1, Dq);   e_or(i) = eps_pred(g2, Dq);
end
sel = [5 9 17 31 33 65 129 257];
fprintf('\n--- Part D: accuracy vs pulse length (M=%d N=%d D=%d) ---\n', Mq, Nq, Dq);
fprintf('  Lg   eps_gs   eps_orth   pulse rate at B=500MHz [G MAC/s]\n');
for L1 = sel
    i = find(Lgv == L1);
    fprintf('%4d  %7.4f  %9.5f   %6.1f\n', L1, e_gs(i), e_or(i), B*L1/1e9);
end
figure('Name','Fig4 accuracy vs cost'); hold on; grid on;
semilogy(Lgv, e_gs, 'o-', 'LineWidth', 1.4); semilogy(Lgv, e_or, 's-', 'LineWidth', 1.4);
set(gca,'YScale','log'); xline(31, '--', 'slide L_g=31', 'Color',[0.7 0.7 0.7], 'LabelOrientation','horizontal', 'HandleVisibility','off');
xlabel('L_g (retained taps)  [cost = B*L_g, 15.5 G MAC/s at 31]'); ylabel('\epsilon  (rms self-interference)');
legend({'Gauss*sinc','orthogonalised'},'Location','northeast');
title('What each tap buys: Nyquist-residual vs pulse length');

% Fig 5: FFT twiddle count: bound vs counted
figure('Name','Fig5 FFT twiddles bound vs counted');
bar(log2(Nlist), [(Nlist/2).*log2(Nlist); nt_all].'); grid on;
xlabel('log_2 N'); ylabel('twiddle multiplications per row');
legend({'slide bound (N/2)log_2N','counted (non-trivial only)'},'Location','northwest');
title('Slide count is conservative');

%% ======================= local functions =======================
function chk(name, val, ref)
    ok = abs(val - ref) <= 1e-3*abs(ref) + 1e-9;
    s = 'FAIL'; if ok, s = 'ok'; end
    fprintf('%-24s %14.6g   slide %14.6g   %s\n', name, val, ref, s);
end

function X = qam16(sz)
    X = ((2*randi([0 3],sz)-3) + 1j*(2*randi([0 3],sz)-3)) / sqrt(10);
end

function y = sincn(t)
    y = ones(size(t));  nz = (t ~= 0);
    y(nz) = sin(pi*t(nz)) ./ (pi*t(nz));
end

function [a, nb, nt] = fft_count(x)      % radix-2 DIT, counts butterflies / nontrivial twiddles
    x = x(:);  N = numel(x);  n = log2(N);
    rev = bin2dec(fliplr(dec2bin(0:N-1, n))) + 1;
    a = x(rev(:));  nb = 0;  nt = 0;
    for st = 1:n
        half = 2^(st-1);
        for b0 = 0:2*half:N-1
            for k = 0:half-1
                w = exp(-1j*pi*k/half);
                i1 = b0 + k + 1;  i2 = b0 + k + half + 1;
                tw = a(i2)*w;  u = a(i1);
                a(i1) = u + tw;  a(i2) = u - tw;
                nb = nb + 1;
                if ~(k == 0 || (st >= 2 && k == half/2)), nt = nt + 1; end
            end
        end
    end
end

function [s, c] = tx_counted(X, g, ell0, D)
    [M, N] = size(X);  P = M*D;  L = N*P;  Lg = numel(g);
    V = zeros(M, N);  c.but = 0;  c.tw = 0;  c.mac = 0;
    for m = 1:M
        [y, nb, nt] = fft_count(conj(X(m,:)).');
        V(m,:) = (conj(y)/N).';              % ifft = conj(fft(conj(x)))/N
        c.but = c.but + nb;  c.tw = c.tw + nt;
    end
    s = zeros(L, 1);
    for q = 0:N-1
        for m = 0:M-1
            v = V(m+1, q+1);
            b0 = D*m + P*q + ell0;
            for i = 1:Lg
                a = mod(b0 + i - 1, L) + 1;
                s(a) = s(a) + v*g(i);
                c.mac = c.mac + 1;
            end
        end
    end
end

function [Y, c] = rx_counted(r, M, N, g, ell0, D)
    P = M*D;  L = N*P;  Lg = numel(g);  r = r(:);
    U = zeros(M, N);  c.mul = 0;  c.add = 0;  c.but = 0;  c.tw = 0;
    for q = 0:N-1
        for m = 0:M-1
            b0 = D*m + P*q + ell0;  acc = 0;
            for i = 1:Lg
                a = mod(b0 + i - 1, L) + 1;
                tm = r(a)*conj(g(i));  c.mul = c.mul + 1;
                if i > 1, acc = acc + tm;  c.add = c.add + 1; else, acc = tm; end
            end
            U(m+1, q+1) = acc;
        end
    end
    Y = zeros(M, N);
    for m = 1:M
        [y, nb, nt] = fft_count(U(m,:).');
        Y(m,:) = y.';  c.but = c.but + nb;  c.tw = c.tw + nt;
    end
end

function ep = eps_pred(g, D)
    Lg = numel(g);
    Rf = conv(g, flipud(conj(g)));
    lags = (-(Lg-1):(Lg-1)).';
    sel = mod(lags, D) == 0;
    kk = lags(sel)/D;
    ep = sqrt(sum(abs(Rf(sel) - double(kk == 0)).^2));
end

function go = orth_full(base, D, L)
    n = numel(base);  c = (n-1)/2;
    full = zeros(L,1);  full(mod(-c:c, L) + 1) = base;
    G = fft(full);  Lq = L/D;  Fo = zeros(Lq,1);
    for k = 0:D-1, Fo = Fo + abs(G(k*Lq + (1:Lq))).^2; end
    go = real(ifft(G .* sqrt(D ./ repmat(Fo, D, 1))));
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