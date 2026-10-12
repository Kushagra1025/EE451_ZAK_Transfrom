%% S5_sampled_vs_continuous.m -- is the sampled receiver (slides p4-6) a good model of the continuous projection?
clear; clc; close all;
rng(5);
M = 16;  N = 8;  alpha = 4;  Wtr = 6;  Tf = M*N;      % all times in units of Ts (Ts = 1)
Dref = 64;                                            % "continuous" reference: Tsamp = Ts/64
X = ((2*randi(4,M,N)-5) + 1j*(2*randi(4,M,N)-5))/sqrt(10);
V = ifft(X, [], 2);                                   % v_m[q]
cn_tr = cnorm(alpha, Wtr, M);                         % continuous normalisation, pulse cut at +-Wtr Ts
cn_un = cnorm(alpha, [],  M);                         % continuous normalisation, untruncated pulse

%% ---- Part A: slides p4, p5, p6 are the same computation ------------------
D = 4;  P = M*D;  L = M*N*D;  Tsamp = 1/D;
r = randn(L,1) + 1j*randn(L,1);                       % any sampled sequence
gsamp = @(ell) gcont(ell/D, alpha, Wtr, M)/cn_tr;     % g[ell] = g(ell*Tsamp), continuous normalisation
U4 = zeros(M, N);                                     % slide p4: Tsamp * sum over ALL samples a
a = (0:L-1).';
for q = 0:N-1
    for m = 0:M-1
        ell = mod(a - D*m - P*q + L/2, L) - L/2;      % a - Dm - Pq, wrapped into [-L/2, L/2)
        U4(m+1,q+1) = Tsamp * sum(r .* conj(gsamp(ell)));
    end
end
Y4 = fft(U4, [], 2);
Lg5  = 2*D*Wtr + 1;  ell5 = (-(Lg5-1)/2 : (Lg5-1)/2).';
gt   = Tsamp * gsamp(ell5);                           % slide p6: g~[ell] = Tsamp * g[ell]
Y56  = zak_sampled(r, M, N, gt, -(Lg5-1)/2, D);       % slides p5/p6: FIR over the Lg retained taps
fprintf('p4 (all samples) vs p5/p6 (FIR, Lg=%d, Tsamp absorbed): max|diff| = %.2e\n', Lg5, max(abs(Y4(:) - Y56(:))));
assert(max(abs(Y4(:) - Y56(:))) < 1e-10);

%% ---- Part B: the Riemann sum, one (m,q) at a time (Figure 1) -------------
m0 = 8;  q0 = 3;  c0 = m0 + q0*M;                     % pulse centre in units of Ts
tf   = c0 + (-8:1/Dref:8).';                          % fine grid around the pulse
integ = @(t) swave(t, V, 0, alpha, Wtr, cn_tr, M, N) .* conj(gcont(wrapt(t - c0, Tf), alpha, Wtr, M)/cn_tr);
If    = integ(tf);
Itrue = sum(If)/Dref;                                 % fine Riemann sum = reference for this (m,q)
Dlist = [1 2 4];
figure('Name','Fig1: sampled projection as a Riemann sum');
for j = 1:3
    Dj = Dlist(j);
    ts = c0 + (-8:1/Dj:8).';
    Is = integ(ts);
    Ssum = sum(Is)/Dj;                                % Tsamp * sum of the samples
    subplot(1,3,j); hold on;
    plot(tf, real(If), 'b-', 'LineWidth', 1.2);
    stem(ts, real(Is), 'r', 'filled', 'MarkerSize', 4);
    grid on; xlabel('t / T_s'); ylabel('Re[ r(t) g^*(t - m_0T_s - q_0T) ]');
    title(sprintf('D=%d: sum = %.4f%+.4fi', Dj, real(Ssum), imag(Ssum)));
    fprintf('D=%d: Tsamp*sum = %.6f%+.6fi    fine reference = %.6f%+.6fi\n', ...
            Dj, real(Ssum), imag(Ssum), real(Itrue), imag(Itrue));
end
sgtitle('blue = integrand (fine grid), red = the samples the ADC-rate sum uses');

%% ---- Part C: convergence with D (Figure 2) --------------------------------
Dv    = [1 2 4 8 16];
cases = {'untruncated, aligned', 'untruncated, offset 0.37 T_s', 'cut at \pm6 T_s, aligned', 'cut at \pm6 T_s, offset 0.37 T_s'};
Wc    = {[], [], Wtr, Wtr};   tauc = [0 0.37 0 0.37];   cnc = [cn_un cn_un cn_tr cn_tr];
errU  = zeros(4, numel(Dv));
for c = 1:4
    Uref = proj_samp(Dref, alpha, Wc{c}, tauc(c), V, cnc(c), M, N);
    for j = 1:numel(Dv)
        Uj = proj_samp(Dv(j), alpha, Wc{c}, tauc(c), V, cnc(c), M, N);
        errU(c,j) = norm(Uj(:) - Uref(:))/norm(Uref(:));
    end
end
disp('relative error of the sampled projection u_m[q] vs the fine reference');
disp('columns: D = 1 2 4 8 16;  rows: untrunc/aligned, untrunc/offset, cut/aligned, cut/offset');
disp(errU);
figure('Name','Fig2: convergence of the sampled projection');
cc = lines(4);  mk = {'o-','s-','o--','s--'};
for c = 1:4
    loglog(Dv, max(errU(c,:), 1e-17), mk{c}, 'Color', cc(c,:), 'LineWidth', 1.4, 'MarkerFaceColor', cc(c,:)); hold on;
end
loglog(Dv(2:end), errU(3,2)*(Dv(2)./Dv(2:end)).^2, 'k:', 'LineWidth', 1.2);
grid on; xticks(Dv); ylim([1e-17 1]);
xlabel('D (samples per T_s)'); ylabel('relative error vs fine reference');
legend([cases, {'1/D^2 guide'}], 'Location', 'southwest');
title('sampled projection vs continuous integral');

Uref3 = proj_samp(Dref, alpha, Wtr, 0, V, cn_tr, M, N);
Dsel = [2 8];
figure('Name','Fig3: where the sampling error sits');
subplot(1,3,1);
imagesc(0:N-1, 0:M-1, abs(Uref3)); axis xy; colorbar;
xlabel('period q'); ylabel('delay m'); title('|u_m[q]| (fine reference)');
for j = 1:2
    Uj = proj_samp(Dsel(j), alpha, Wtr, 0, V, cn_tr, M, N);
    subplot(1,3,1+j);
    imagesc(0:N-1, 0:M-1, 20*log10(abs(Uj - Uref3)/max(abs(Uref3(:))) + 1e-17)); axis xy; colorbar;
    xlabel('period q'); ylabel('delay m'); title(sprintf('error in dB re max|u|, D=%d (cut pulse)', Dsel(j)));
end

%% ---- Part D: what S4's discrete model means in continuous time (Figure 4) -
epsC_tr = eps_cont(alpha, Wtr, M);   epsC_un = eps_cont(alpha, [], M);
fprintf('continuous-time Nyquist residual eps: untruncated = %.4f, cut at +-%d Ts = %.4f\n', epsC_un, Wtr, epsC_tr);
Dd = [1 2 4 8 16];  epsD = zeros(size(Dd));
for j = 1:numel(Dd)
    [gj, ~] = pulse_gs(Dd(j), 2*Dd(j)*Wtr + 1, alpha, M);
    [~,~,~,~,epsD(j)] = pulse_stats(gj, Dd(j));
    fprintf('D=%2d: discrete (S4) eps_pred = %.5f\n', Dd(j), epsD(j));
end
figure('Name','Fig4: discrete model vs continuous time');
plot(Dd, epsD, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto'); hold on;
yline(epsC_tr, 'r--', 'LineWidth', 1.3);
yline(epsC_un, 'g--', 'LineWidth', 1.3);
set(gca, 'XScale', 'log'); xticks(Dd); grid on;
xlabel('D'); ylabel('\epsilon (Nyquist residual)');
legend('discrete model (S4), cut at \pm6 T_s', 'continuous time, cut at \pm6 T_s', 'continuous time, untruncated', 'Location', 'southeast');
title('the discrete model converges to the continuous pulse, except at D = 1');

% why D = 1 is a trap: continuous normalisation vs sampled normalisation
fprintf('integral of |g|^2 dt for the untruncated pulse (before normalising) = %.4f\n', cn_un^2);
Ud1 = proj_samp(1, alpha, [], 0, V, cn_un, M, N);
Yd1 = fft(Ud1, [], 2);
fprintf('D=1 sampled chain, continuous normalisation: eps = %.4f   (1/cn^2 - 1 = %.4f)\n', ...
        norm(Yd1(:) - X(:))/norm(X(:)), 1/cn_un^2 - 1);

%% ==================== local functions (must stay at the end) ==============
function y = sincn(x)
    y = ones(size(x));  nz = x ~= 0;
    y(nz) = sin(pi*x(nz))./(pi*x(nz));
end

function g = gcont(t, alpha, W, M)
% continuous Gaussian x sinc (slide p15); tau_s = Ts, T = M Ts; W = [] means untruncated
    g = exp(-pi*alpha*(t/M).^2) .* sincn(t);
    if ~isempty(W), g = g .* (abs(t) <= W); end
end

function d = wrapt(t, Tf)
    d = mod(t + Tf/2, Tf) - Tf/2;                       % circular frame of length Tf
end

function c = cnorm(alpha, W, M)
% sqrt( integral |g|^2 dt ) on a fine grid
    dt = 1/256;  t = (-40:dt:40).';
    g = gcont(t, alpha, W, M);
    c = sqrt(sum(g.^2)*dt);
end

function e = eps_cont(alpha, W, M)
% Nyquist residual of the continuous pulse: sqrt( sum_k |R(k Ts) - delta_k|^2 ), R = autocorrelation integral
    dt = 1/256;  t = (-60:dt:60-dt).';
    g = gcont(t, alpha, W, M);  g = g/sqrt(sum(g.^2)*dt);
    e2 = 0;
    for k = -60:60
        sh = k/dt;                                      % integer number of grid points
        R  = sum(g(1+max(0,sh):end+min(0,sh)) .* g(1+max(0,-sh):end-max(0,sh)))*dt;
        e2 = e2 + abs(R - (k == 0))^2;
    end
    e = sqrt(e2);
end

function s = swave(t, V, tau0, alpha, W, cn, M, N)
% continuous transmit waveform s(t - tau0) = sum_{m,q} v_m[q] g(t - tau0 - (m + qM)Ts)/cn, evaluated analytically
    Tf = M*N;  cen = (0:M-1).' + M*(0:N-1);
    s = zeros(size(t));
    for i = 1:M*N
        s = s + V(i) * gcont(wrapt(t - tau0 - cen(i), Tf), alpha, W, M) / cn;
    end
end

function U = proj_samp(D, alpha, W, tau0, V, cn, M, N)
% sampled projection u_m[q] = Tsamp * sum_a r[a] g*[a - Dm - Pq]  with r(t) = s(t - tau0) sampled every Ts/D
    Tf = M*N;  L = M*N*D;  dt = 1/D;
    t   = (0:L-1).' * dt;
    cen = (0:M-1).' + M*(0:N-1);
    s = swave(t, V, tau0, alpha, W, cn, M, N);
    U = zeros(M, N);
    for i = 1:M*N
        gg = gcont(wrapt(t - cen(i), Tf), alpha, W, M) / cn;
        U(i) = sum(s .* conj(gg)) * dt;
    end
end

function [g, ell0] = pulse_gs(D, Lg, alpha, M)
    c   = (Lg-1)/2;
    ell = (-c:c).';
    t   = ell/D;
    g   = exp(-pi*alpha*(t/M).^2) .* sincn(t);
    g   = g / sqrt(sum(abs(g).^2));
    ell0 = -c;
end

function [lags, Rfull, kk, Rk, ep] = pulse_stats(g, D)
    Lg    = numel(g);
    Rfull = conv(g, flipud(conj(g)));
    lags  = (-(Lg-1):(Lg-1)).';
    sel   = mod(lags, D) == 0;
    kk    = lags(sel)/D;
    Rk    = Rfull(sel);
    ep    = sqrt(sum(abs(Rk - (kk == 0)).^2));
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