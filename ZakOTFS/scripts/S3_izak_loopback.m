%% S3_izak_loopback.m  -- inverse Zak transform and Level-0 loopback (slides p11-13, 17, 21)
clear; clc; close all;
rng(3);
M = 16;  N = 8;

%% ---- Part A: step 1 of the IZak = IFFT over n (Figure 1) ----------------
% One DD point X[m0,n0] = 1  ->  v_m0[q] = (1/N) exp(+j 2 pi n0 q / N)
m0 = 3;  q = 0:N-1;
cols = lines(4);
figure('Name','Fig1: Doppler phase trajectories');
subplot(1,3,1); hold on;
th = linspace(0, 2*pi, 200);
plot(cos(th)/N, sin(th)/N, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
for n0 = 0:3
    X0 = zeros(M, N);  X0(m0+1, n0+1) = 1;
    V0 = ifft(X0, [], 2);
    vth = exp(1j*2*pi*n0*q/N)/N;
    assert(max(abs(V0(m0+1,:) - vth)) < 1e-14);
    plot(real(V0(m0+1,:)), imag(V0(m0+1,:)), 'o', 'Color', cols(n0+1,:), ...
         'MarkerFaceColor', cols(n0+1,:), 'MarkerSize', 14 - 3*n0, 'DisplayName', sprintf('n_0=%d', n0));
end
axis equal; grid on; xlabel('Re v'); ylabel('Im v'); title('v_{m_0}[q] in the complex plane'); legend('show');
subplot(1,3,2); hold on;
for n0 = 0:3
    vth = exp(1j*2*pi*n0*q/N)/N;
    plot(q, real(vth), 'o-', 'Color', cols(n0+1,:), 'LineWidth', 1.2, 'DisplayName', sprintf('n_0=%d', n0));
end
grid on; xlabel('period index q'); ylabel('Re v_{m_0}[q]'); title('real part vs q'); legend('show');
subplot(1,3,3);
Xrow = zeros(M, N);  Xrow(m0+1, 2) = 1;  Xrow(m0+1, 4) = 0.7;       % Doppler bins n=1 and n=3
Vrow = ifft(Xrow, [], 2);
plot(q, real(Vrow(m0+1,:)), 'ko-', 'LineWidth', 1.5); hold on;
plot(q, real(exp(1j*2*pi*1*q/N)/N + 0.7*exp(1j*2*pi*3*q/N)/N), 'r.', 'MarkerSize', 16);
grid on; xlabel('period q'); ylabel('Re v_{m_0}[q]');
title('two Doppler bins: IFFT (black) = sum of trajectories (red)'); legend('IFFT', 'sum of two exponentials');

%% ---- Part B: Level-0 loopback, D = 1, one-sample pulse (Figure 2) -------
X = ((2*randi(4,M,N)-5) + 1j*(2*randi(4,M,N)-5))/sqrt(10);   % random 16-QAM-like grid, unit average power
D = 1;  g1 = 1;  ell0 = 0;
s1 = izak_sampled(X, g1, ell0, D);
Y1 = zak_sampled(s1, M, N, g1, ell0, D);
err1 = max(abs(Y1(:) - X(:)));
fprintf('Level 0 (D=1): max|Y - X| = %.2e\n', err1);
assert(err1 < 1e-12);
assert(max(abs(s1 - reshape(ifft(X,[],2), [], 1))) < 1e-12);   % for D=1, s is just V(:)

figure('Name','Fig2: Level-0 loopback (D=1)');
subplot(2,2,1);
plot(real(X(:)), imag(X(:)), 'o'); axis equal; grid on;
xlabel('Re'); ylabel('Im'); title('transmitted X[m,n]');
subplot(2,2,2); hold on;
cq = lines(N);
for qq = 0:N-1
    kk = (0:M-1) + qq*M;
    stem(kk, real(s1(kk+1)), 'Color', cq(qq+1,:), 'MarkerFaceColor', cq(qq+1,:), 'MarkerSize', 3);
end
xlim([-1 M*N]); xlabel('sample index a'); ylabel('Re s[a]'); title('IZak output (colour = period q)');
subplot(2,2,3);
plot(real(X(:)), imag(X(:)), 'o', 'MarkerSize', 10); hold on;
plot(real(Y1(:)), imag(Y1(:)), 'r.', 'MarkerSize', 12);
axis equal; grid on; xlabel('Re'); ylabel('Im'); title('recovered Y (red) on top of X (circles)');
subplot(2,2,4);
imagesc(0:N-1, 0:M-1, abs(Y1 - X)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|Y - X| per bin');

%% ---- Part C: rectangular pulse of width Ts = sample-and-hold (Figure 3) -
D = 4;  P = M*D;
g = ones(D,1)/sqrt(D);  ell0 = 0;            % unit-energy rectangular pulse, D samples wide
s4 = izak_sampled(X, g, ell0, D);
Y4 = zak_sampled(s4, M, N, g, ell0, D);
fprintf('Rect pulse, D=4: max|Y - X| = %.2e\n', max(abs(Y4(:) - X(:))));
assert(max(abs(s4*sqrt(D) - repelem(s1, D))) < 1e-12);    % each D=1 sample is held for D ticks

figure('Name','Fig3: rectangular pulse = sample-and-hold');
subplot(2,1,1); hold on;
stairs(0:2*P-1, real(s4(1:2*P))*sqrt(D), 'LineWidth', 1.2);
stem(D*(0:2*M-1), real(s1(1:2*M)), 'r', 'filled', 'MarkerSize', 4);
xline(P, 'k--'); xlabel('DAC sample index a'); ylabel('Re');
title('first two Zak periods: s_{D=4} \cdot \surd D (blue) and the D=1 samples (red)');
subplot(2,1,2); hold on;
stairs(0:2*P-1, imag(s4(1:2*P))*sqrt(D), 'LineWidth', 1.2);
stem(D*(0:2*M-1), imag(s1(1:2*M)), 'r', 'filled', 'MarkerSize', 4);
xline(P, 'k--'); xlabel('DAC sample index a'); ylabel('Im');

%% ---- Part D: loopback error and energy versus D (Figure 4) --------------
Dv = [1 2 4 8 16];
errD = zeros(size(Dv));  ratioE = zeros(size(Dv));
for i = 1:numel(Dv)
    Di = Dv(i);  gi = ones(Di,1)/sqrt(Di);
    si = izak_sampled(X, gi, 0, Di);
    Yi = zak_sampled(si, M, N, gi, 0, Di);
    errD(i)   = max(abs(Yi(:) - X(:)));
    ratioE(i) = sum(abs(si).^2) / (sum(abs(X(:)).^2)/N);     % IZak shrinks energy by N
    fprintf('D=%2d: max|Y-X| = %.2e,  E_s/(E_X/N) = %.6f\n', Di, errD(i), ratioE(i));
end
figure('Name','Fig4: loopback error and energy vs D');
subplot(1,2,1);
semilogy(Dv, max(errD, 1e-17), 'o-', 'LineWidth', 1.5); set(gca, 'XScale', 'log'); grid on;
xticks(Dv); ylim([1e-18 1e-12]);
xlabel('D (DAC samples per delay bin)'); ylabel('max|Y - X|'); title('loopback error, rectangular pulse');
subplot(1,2,2);
plot(Dv, ratioE, 'o-', 'LineWidth', 1.5); set(gca, 'XScale', 'log'); grid on;
xticks(Dv); ylim([0.9 1.1]);
xlabel('D'); ylabel('E_s / (E_X / N)'); title('energy check (should be 1)');

%% ---- Part E: three IZak implementations must agree (Figure 5) -----------
D = 4;  P = M*D;  L = N*P;  g = ones(D,1)/sqrt(D);
sA = izak_sampled(X, g, 0, D);                         % (a) fast: IFFT + pulse placement
Fm   = exp(1j*2*pi*(0:N-1).'*(0:N-1)/N);               % Fm(n+1,q+1) = e^{+j 2 pi n q / N}
Vexp = X*Fm/N;                                         % explicit N-term sum instead of ifft
sB = zeros(L,1);                                       % (b) literal slide formula
for m = 0:M-1
    for qq = 0:N-1
        for i = 1:D
            a = mod(D*m + P*qq + i - 1, L);
            sB(a+1) = sB(a+1) + Vexp(m+1,qq+1)*g(i);
        end
    end
end
sC = repelem(reshape(ifft(X,[],2), [], 1), D)/sqrt(D); % (c) Level-0 signal, held D ticks
eAB = max(abs(sA - sB));  eAC = max(abs(sA - sC));
fprintf('IZak: |fast - literal| = %.2e,  |fast - hold| = %.2e\n', eAB, eAC);
assert(eAB < 1e-12 && eAC < 1e-12);
figure('Name','Fig5: three IZak implementations');
subplot(1,2,1);
semilogy(1:2, max([eAB eAC], 1e-18), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); grid on;
xlim([0.5 2.5]); ylim([1e-19 1e-12]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'fast vs literal', 'fast vs hold'});
ylabel('max|difference|'); title('same signal, three ways');
subplot(1,2,2);
imagesc(0:N-1, 0:P-1, abs(reshape(sA - sB, P, N))); axis xy; colorbar;
xlabel('period q'); ylabel('sample within period'); title('|fast - literal| per sample');

%% ---- Part F: loopback over many random frames (Figure 6) ----------------
nFrames = 300;  Dl = [1 4 8];
E = zeros(nFrames, numel(Dl));
for f = 1:nFrames
    Xf = ((2*randi(4,M,N)-5) + 1j*(2*randi(4,M,N)-5))/sqrt(10);
    for i = 1:numel(Dl)
        gi = ones(Dl(i),1)/sqrt(Dl(i));
        sf = izak_sampled(Xf, gi, 0, Dl(i));
        E(f,i) = max(max(abs(zak_sampled(sf, M, N, gi, 0, Dl(i)) - Xf)));
    end
end
fprintf('%d frames x %d values of D: worst loopback error = %.2e\n', nFrames, numel(Dl), max(E(:)));
figure('Name','Fig6: random-frame loopback errors');
hold on;
for i = 1:numel(Dl)
    histogram(log10(max(E(:,i), 1e-20)), 'BinWidth', 0.1, 'DisplayName', sprintf('D=%d', Dl(i)));
end
grid on; legend('show'); xlabel('log_{10} max|Y - X|'); ylabel('number of frames');
title(sprintf('%d random frames', nFrames));

%% ---- Part G: a sparse DD grid becomes a waveform and comes back (Fig 7) -
D = 4;  P = M*D;  g = ones(D,1)/sqrt(D);
Xs = zeros(M, N);
Xs(3+1, 1+1) = 1;  Xs(3+1, 5+1) = 1j;  Xs(10+1, 2+1) = -1;  Xs(13+1, 6+1) = 0.7;
ss = izak_sampled(Xs, g, 0, D);
Ys = zak_sampled(ss, M, N, g, 0, D);
fprintf('Sparse grid: max|Y - X| = %.2e\n', max(abs(Ys(:) - Xs(:))));
assert(max(abs(Ys(:) - Xs(:))) < 1e-12);
figure('Name','Fig7: sparse DD grid -> waveform -> DD grid');
subplot(1,3,1);
imagesc(0:N-1, 0:M-1, abs(Xs)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|X[m,n]| (4 occupied cells)');
subplot(1,3,2);
imagesc(0:N-1, 0:P-1, abs(reshape(ss, P, N))); axis xy; colorbar;
xlabel('period q'); ylabel('sample within period (a mod P)'); title('|s[a]| folded into periods');
subplot(1,3,3);
imagesc(0:N-1, 0:M-1, abs(Ys)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|Y[m,n]| after Zak');

%% ==================== local functions (must stay at the end) ==============
function s = izak_sampled(X, g, ell0, D)
% s[a] = sum_m sum_q v_m[q] g[a - D m - P q],  v_m[q] = (1/N) sum_n X[m,n] e^{+j2pi nq/N}
% g(i) is the pulse tap at relative index ell = ell0 + i - 1.  Frame is treated as circular.
    [M, N] = size(X);
    P = M*D;  L = N*P;
    V = ifft(X, [], 2);                          % step 1: IFFT over n (includes the 1/N)
    up = zeros(L, 1);
    k = (0:M-1).' + M*(0:N-1);                   % k(m+1,q+1) = m + q*M
    up(D*k(:) + 1) = V(:);                       % impulses at a = D m + P q = D (m + qM)
    s = zeros(L, 1);
    for i = 1:numel(g)                           % step 2: pulse placement
        s = s + g(i) * circshift(up, ell0 + i - 1);
    end
end

function Y = zak_sampled(r, M, N, g, ell0, D)
% u_m[q] = sum_ell r[D m + P q + ell] conj(g[ell]);  Y[m,n] = sum_q u_m[q] e^{-j2pi nq/N}
    L = M*N*D;  r = r(:);
    acc = zeros(L, 1);
    for i = 1:numel(g)
        acc = acc + circshift(r, -(ell0 + i - 1)) * conj(g(i));    % r[a0+ell] * conj(g[ell])
    end
    k = (0:M-1).' + M*(0:N-1);
    U = reshape(acc(D*k(:) + 1), M, N);          % u_m[q]
    Y = fft(U, [], 2);                           % N-point FFT over q
end