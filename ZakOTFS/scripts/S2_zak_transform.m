%% S2_zak_transform.m  -- Zak transform = reshape + FFT (slides p2-3)
clear; clc; close all;
B = 500e6;                      % only used to print physical units

%% ---- Part A: fold and FFT, step by step (Figure 1) ----------------------
rng(0);
M = 16;  N = 8;
r = randn(M*N,1) + 1j*randn(M*N,1);      % random received samples, r(k+1) = r[k]

R = reshape(r, M, N);                    % fold: column q+1 = period q, row m+1 = delay m
Y = fft(R, [], 2);                       % FFT across q

figure('Name','Fig1: Reshape + FFT');
subplot(1,3,1);
stem(0:M*N-1, real(r), 'filled', 'MarkerSize', 3); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('Re r[k]');
title('1. samples, cut every M');
subplot(1,3,2);
imagesc(0:N-1, 0:M-1, real(R)); axis xy; colorbar;
xlabel('period q'); ylabel('delay m'); title('2. R[m,q] = r[m+qM]');
subplot(1,3,3);
imagesc(0:N-1, 0:M-1, abs(Y)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('3. |Y[m,n]|');

%% ---- Part B: three implementations must agree (Figure 2) ----------------
Mb = 64;  Nb = 64;
rb = randn(Mb*Nb,1) + 1j*randn(Mb*Nb,1);

Y1 = dzt_bruteforce(rb, Mb, Nb);                       % slide formula, triple loop
Y2 = dzt(rb, Mb, Nb);                                  % reshape + FFT
W  = exp(-1j*2*pi*(0:Nb-1).'*(0:Nb-1)/Nb);             % DFT matrix, W(q+1,n+1) = e^{-j2pi nq/N}
Y3 = reshape(rb, Mb, Nb) * W;                          % FFT as a matrix product

e12 = max(abs(Y1(:) - Y2(:)));
e32 = max(abs(Y3(:) - Y2(:)));
fprintf('max|brute - fft|  = %.2e\n', e12);
fprintf('max|matrix - fft| = %.2e\n', e32);
fprintf('sum|Y|^2 / sum|r|^2 = %.4f   (N = %d)\n', sum(abs(Y2(:)).^2)/sum(abs(rb).^2), Nb);
assert(e12 < 1e-9 && e32 < 1e-9);

figure('Name','Fig2: Three implementations');
tiledlayout(1,2);
nexttile;
semilogy(1:2, [e12 e32], 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); grid on;
xlim([0.5 2.5]); set(gca, 'XTick', 1:2, 'XTickLabel', {'brute vs FFT','matrix vs FFT'});
ylabel('max abs difference'); title('same transform, three ways');
nexttile;
imagesc(abs(Y1 - Y2)); axis xy; colorbar;
xlabel('n'); ylabel('m'); title('|brute - FFT| per bin');

%% ---- Part C: an impulse gives a row (Figure 3) --------------------------
M = 16;  N = 8;  m0 = 5;  q0 = 2;
r = zeros(M*N,1);
r(m0 + q0*M + 1) = 1;
Y = dzt(r, M, N);

Ytheory = zeros(M, N);
Ytheory(m0+1, :) = exp(-1j*2*pi*(0:N-1)*q0/N);
fprintf('Impulse: max|Y - theory| = %.2e\n', max(abs(Y(:) - Ytheory(:))));
assert(max(abs(Y(:) - Ytheory(:))) < 1e-12);

figure('Name','Fig3: Impulse -> row');
subplot(1,3,1);
stem(0:M*N-1, abs(r), 'filled'); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('|r[k]|');
title(sprintf('impulse at m_0=%d, q_0=%d', m0, q0));
subplot(1,3,2);
imagesc(0:N-1, 0:M-1, abs(Y)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|Y[m,n]|: one row');
subplot(1,3,3);
stem(0:N-1, angle(Y(m0+1,:)), 'filled'); hold on;
plot(0:N-1, angle(Ytheory(m0+1,:)), 'ro', 'MarkerSize', 10);
ylim([-pi pi]); grid on; xlabel('Doppler bin n'); ylabel('phase (rad)');
title('phase along the row (red = theory)');

%% ---- Part D: a tone gives a column (Figure 4) ---------------------------
n0 = 3;  k = (0:M*N-1).';
r = exp(1j*2*pi*n0*k/(M*N));
Y = dzt(r, M, N);

Ytheory = zeros(M, N);
Ytheory(:, n0+1) = N*exp(1j*2*pi*n0*(0:M-1).'/(M*N));
fprintf('Tone: max|Y - theory| = %.2e\n', max(abs(Y(:) - Ytheory(:))));
fprintf('Tone offset = n0*B/(M*N) = %.5f MHz\n', n0*B/(M*N)/1e6);
assert(max(abs(Y(:) - Ytheory(:))) < 1e-10);

figure('Name','Fig4: Tone -> column');
subplot(1,3,1);
stem(0:M*N-1, real(r), 'filled', 'MarkerSize', 3); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('Re r[k]');
title(sprintf('tone, n_0 = %d', n0));
subplot(1,3,2);
imagesc(0:N-1, 0:M-1, abs(Y)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|Y[m,n]|: one column');
subplot(1,3,3);
plot(0:M-1, angle(Y(:,n0+1)), 'o', 'MarkerSize', 8); hold on;
plot(0:M-1, angle(Ytheory(:,n0+1)), 'r.', 'MarkerSize', 12);
grid on; xlabel('delay m'); ylabel('phase (rad)'); title('phase down the column (red = theory)');

%% ---- Part E: one DD point = phase-rotating pulse train (Figure 5) -------
m0 = 11;  n0 = 2;
r = zeros(M*N,1);
for q = 0:N-1
    r(m0 + q*M + 1) = exp(1j*2*pi*n0*q/N);     % one pulse per period, phase rotating with q
end
Y = dzt(r, M, N);

Ytheory = zeros(M, N);  Ytheory(m0+1, n0+1) = N;
fprintf('DD point: max|Y - theory| = %.2e\n', max(abs(Y(:) - Ytheory(:))));
assert(max(abs(Y(:) - Ytheory(:))) < 1e-10);

figure('Name','Fig5: Pulse train -> DD point');
subplot(1,3,1);
stem(0:M*N-1, real(r), 'filled'); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('Re r[k]'); title('real part of the pulse train');
subplot(1,3,2);
stem(0:M*N-1, imag(r), 'filled'); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('Im r[k]'); title('imaginary part');
subplot(1,3,3);
imagesc(0:N-1, 0:M-1, abs(Y)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title(sprintf('|Y|: one point (m=%d, n=%d)', m0, n0));

%% ---- Part F: superposition with noise (Figure 6) ------------------------
rng(2);  sigma = 0.5;
k = (0:M*N-1).';
r = zeros(M*N,1);
r(3 + 1*M + 1) = N;                                   % impulse -> row m=3
r = r + exp(1j*2*pi*5*k/(M*N));                       % tone -> column n=5
for q = 0:N-1
    r(11 + q*M + 1) = r(11 + q*M + 1) + exp(1j*2*pi*2*q/N);   % DD point (11,2)
end
r = r + sigma*(randn(M*N,1) + 1j*randn(M*N,1))/sqrt(2);       % complex noise, variance sigma^2
Y = dzt(r, M, N);

mask = true(M, N);  mask(3+1,:) = false;  mask(:,5+1) = false;  mask(11+1,2+1) = false;
fprintf('Noise-only bins: mean|Y|^2 = %.2f   (expected N*sigma^2 = %.2f)\n', ...
        mean(abs(Y(mask)).^2), N*sigma^2);

figure('Name','Fig6: Superposition');
subplot(1,2,1);
stem(0:M*N-1, real(r), 'filled', 'MarkerSize', 3); hold on;
for q = 0:N, xline(q*M - 0.5, 'r--'); end
xlim([-1 M*N]); xlabel('sample index k'); ylabel('Re r[k]'); title('received signal (looks like noise)');
subplot(1,2,2);
imagesc(0:N-1, 0:M-1, abs(Y)); axis xy; colorbar;
xlabel('Doppler bin n'); ylabel('delay m'); title('|Y[m,n]|: row + column + point + noise');

%% ---- Part G: why the FFT matters (Figure 7) -----------------------------
Mt = 64;  Nv = 2.^(4:11);
tF = zeros(size(Nv));  tD = zeros(size(Nv));
for i = 1:numel(Nv)
    Nt = Nv(i);
    Rt = randn(Mt, Nt) + 1j*randn(Mt, Nt);
    Wt = exp(-1j*2*pi*(0:Nt-1).'*(0:Nt-1)/Nt);
    tF(i) = timeit(@() fft(Rt, [], 2));
    tD(i) = timeit(@() Rt*Wt);
end
cntD = Mt*Nv.^2;               % direct: M N^2 complex MACs (slide p17)
cntF = Mt*(Nv/2).*log2(Nv);    % FFT:    M (N/2) log2 N butterflies (slide p17)

figure('Name','Fig7: FFT vs direct');
tiledlayout(1,2);
nexttile;
loglog(Nv, cntD, 'o-', Nv, cntF, 's-', 'LineWidth', 1.5); grid on;
legend('direct  M N^2', 'FFT  M (N/2) log_2 N', 'Location', 'northwest');
xlabel('N'); ylabel('operations'); title('operation counts (formulas)');
nexttile;
loglog(Nv, tD*1e3, 'o-', Nv, tF*1e3, 's-', 'LineWidth', 1.5); grid on;
legend('matrix product (direct)', 'fft', 'Location', 'northwest');
xlabel('N'); ylabel('time (ms)'); title(sprintf('measured time, M = %d', Mt));
fprintf('Formula ratio direct/FFT at N=128: %.1f\n', (Mt*128^2)/(Mt*64*log2(128)));

%% ==================== local functions (must stay at the end) ==============
function Y = dzt(r, M, N)
% Y(m+1,n+1) = sum_q r(m+q*M+1) * exp(-j*2*pi*n*q/N)   (D = 1, impulse receive pulse)
    R = reshape(r(:), M, N);     % column q+1 = period q, row m+1 = delay m
    Y = fft(R, [], 2);           % N-point FFT across q, for every delay row
end

function Y = dzt_bruteforce(r, M, N)
% The slide formula, literally, with three nested loops.
    Y = zeros(M, N);
    for m = 0:M-1
        for n = 0:N-1
            acc = 0;
            for q = 0:N-1
                acc = acc + r(m + q*M + 1) * exp(-1j*2*pi*n*q/N);
            end
            Y(m+1, n+1) = acc;
        end
    end
end