%% S4_pulses.m -- pulse shaping: normalisation, Nyquist condition, ISI operator H (slides p3, p5-6, p15)
clear; clc; close all;
rng(4);
M = 16;  N = 8;  D = 4;  L = M*N*D;
X = ((2*randi(4,M,N)-5) + 1j*(2*randi(4,M,N)-5))/sqrt(10);      % random 16-QAM-like grid

%% ---- Part A: the pulses (Figure 1) --------------------------------------
Lg = 65;                                          % 65 taps = +-8 Ts at D = 4
[gR,   e0R  ] = pulse_rect(D);
[gG4,  e0G4 ] = pulse_gs(D, Lg, 4,  M);           % Gaussian x sinc, alpha = 4
[gG16, e0G16] = pulse_gs(D, Lg, 16, M);           % Gaussian x sinc, alpha = 16
[gO,   e0O  ] = pulse_orth(D, Lg, 16, M, L);      % orthogonalised version of the alpha = 16 pulse
names = {'rect (width T_s)', 'Gauss x sinc, \alpha=4', 'Gauss x sinc, \alpha=16', 'orthogonalised, \alpha=16'};
G  = {gR, gG4, gG16, gO};
E0 = [e0R, e0G4, e0G16, e0O];
for i = 1:4
    assert(abs(sum(abs(G{i}).^2) - 1) < 1e-12);   % slide p3: sum |g[a]|^2 = 1
end
figure('Name','Fig1: pulses');
for i = 1:4
    subplot(2,2,i); hold on;
    ell = E0(i) + (0:numel(G{i})-1).';
    t = ell/D;                                    % time in units of Ts
    plot(t, G{i}, '-', 'LineWidth', 1.2);
    isTs = mod(ell, D) == 0;
    plot(t(isTs), G{i}(isTs), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    xline(0, 'k:'); grid on; xlim([-8.5 8.5]);
    xlabel('t / T_s'); ylabel('g[\ell]'); title(names{i});
end
sgtitle('red = samples at multiples of T_s (what the Nyquist test looks at)');

%% ---- Part B: Nyquist test on the pulse autocorrelation (Figure 2) --------
figure('Name','Fig2: autocorrelation at multiples of Ts');
for i = 1:4
    [lags, Rfull, kk, Rk, ep] = pulse_stats(G{i}, D);
    off = Rk(kk ~= 0);
    fprintf('%-28s R[0]=%.4f  max|R[kD]|(k~=0)=%.2e  eps_pred=%.4f\n', ...
            names{i}, Rk(kk==0), max([abs(off); 0]), ep);
    subplot(2,2,i); hold on;
    plot(lags/D, Rfull, '-', 'LineWidth', 1.2);
    plot(kk, Rk, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    grid on; xlabel('lag / T_s'); ylabel('R[\ell]');
    title(sprintf('%s,  \\epsilon = %.4f', names{i}, ep));
end
sgtitle('red circles = R[kD]: the Nyquist condition wants 1 at k=0 and 0 elsewhere');

%% ---- Part C: the end-to-end operator H (Figure 3) -----------------------
[gA, eA] = pulse_rect(D);
[gB, eB] = pulse_gs(D, 33, 4, M);
[gC, eC] = pulse_orth(D, 65, 16, M, L);
Hn = {'rect', 'Gauss x sinc, \alpha=4, 33 taps', 'orthogonalised, \alpha=16, 65 taps'};
Gs = {gA, gB, gC};  Es = [eA eB eC];
Hs = cell(1,3);
figure('Name','Fig3: end-to-end operator H');
for i = 1:3
    Hs{i} = build_H(Gs{i}, Es(i), D, M, N);
    [~,~,~,~,ep] = pulse_stats(Gs{i}, D);
    fro = norm(Hs{i} - eye(M*N), 'fro')/sqrt(M*N);
    fprintf('%-34s ||H-I||_F/sqrt(MN) = %.4f   eps_pred = %.4f\n', Hn{i}, fro, ep);
    assert(abs(fro - ep) < 1e-9);
    subplot(2,3,i);
    imagesc(abs(Hs{i})); axis image; colorbar; caxis([0 1]);
    xlabel('input index m + M n'); ylabel('output index'); title(['|H|: ' Hn{i}]);
    subplot(2,3,3+i);
    imagesc(20*log10(abs(Hs{i} - eye(M*N)) + 1e-16)); axis image; colorbar; caxis([-80 0]);
    xlabel('input index m + M n'); ylabel('output index'); title('|H - I| in dB');
end

%% ---- Part D: structure of H (Figure 4) ----------------------------------
H = Hs{2};  [~, ~, kk, Rk] = pulse_stats(gB, D);
m0 = 8;  n0 = 3;
col = reshape(H(:, m0 + M*n0 + 1), M, N);            % response to a unit symbol at (m0, n0)
figure('Name','Fig4: structure of H');
subplot(1,3,1);
imagesc(abs(H(1:M, 1:M))); axis image; colorbar;
xlabel('input delay m'); ylabel('output delay m'); title('block n = 0 of |H| (mixing along delay)');
subplot(1,3,2);
stem((0:M-1) - m0, real(col(:, n0+1)), 'filled'); hold on;
sel = abs(kk) <= 7;
plot(kk(sel), Rk(sel), 'ro', 'MarkerSize', 9);
grid on; xlabel('output delay m - m_0'); ylabel('H'); legend('H column', 'autocorrelation R[kD]');
title('response to one symbol');
subplot(1,3,3);
energyPerN = sum(abs(col).^2, 1);                    % energy landing in each Doppler bin
semilogy(0:N-1, max(energyPerN, 1e-40), 'o-', 'LineWidth', 1.5); grid on;
xlabel('output Doppler bin n'); ylabel('energy'); title(sprintf('leakage stays in n = %d', n0));
fprintf('Doppler leakage (energy outside n=%d): %.2e   (energy inside: %.4f)\n', ...
        n0, sum(energyPerN([1:n0, n0+2:N])), energyPerN(n0+1));

%% ---- Part E: loopback error versus pulse span, alpha, D (Figure 5) ------
Wv  = [1 2 3 4 6 8];                                 % half-span in units of Ts
T   = 30;                                            % random frames per point
fam = {'plain sinc (\alpha=0)', 'Gauss x sinc \alpha=4', 'Gauss x sinc \alpha=16', 'orthogonalised \alpha=16'};
Eprd = zeros(numel(fam), numel(Wv));  Emea = Eprd;
for f = 1:numel(fam)
    for j = 1:numel(Wv)
        Lgj = 2*D*Wv(j) + 1;
        switch f
            case 1, [g, e0] = pulse_gs(D, Lgj, 0,  M);
            case 2, [g, e0] = pulse_gs(D, Lgj, 4,  M);
            case 3, [g, e0] = pulse_gs(D, Lgj, 16, M);
            case 4, [g, e0] = pulse_orth(D, Lgj, 16, M, L);
        end
        [~,~,~,~,Eprd(f,j)] = pulse_stats(g, D);
        Emea(f,j) = meas_eps(g, e0, D, M, N, T);
    end
end
disp('half-span W (Ts):'); disp(Wv);
disp('predicted eps (rows: plain sinc, a=4, a=16, orthogonalised):'); disp(Eprd);
disp('measured  eps:'); disp(Emea);

Dl = [1 2 4 8];  Pd = zeros(size(Dl));  Md = Pd;     % effect of oversampling D at half-span 4 Ts
for j = 1:numel(Dl)
    [g, e0] = pulse_gs(Dl(j), 2*Dl(j)*4 + 1, 4, M);
    [~,~,~,~,Pd(j)] = pulse_stats(g, Dl(j));
    Md(j) = meas_eps(g, e0, Dl(j), M, N, T);
end
fprintf('D sweep (Gauss x sinc, alpha=4, half-span 4 Ts): pred = %s\n', mat2str(Pd, 4));
fprintf('                                                 meas = %s\n', mat2str(Md, 4));

figure('Name','Fig5: loopback error vs span, alpha, D');
cc = lines(numel(fam));
subplot(1,3,1); hold on;
for f = 1:numel(fam)
    semilogy(Wv, Eprd(f,:), '-', 'Color', cc(f,:), 'LineWidth', 1.3, 'DisplayName', fam{f});
    semilogy(Wv, Emea(f,:), 'o', 'Color', cc(f,:), 'MarkerFaceColor', cc(f,:), 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log'); grid on; legend('show', 'Location', 'southwest');
xlabel('pulse half-span (T_s)'); ylabel('\epsilon = rms(Y-X)/rms(X)'); title('lines = prediction, dots = measured');
subplot(1,3,2);
semilogy(Dl, max(Pd, 1e-17), '-', 'LineWidth', 1.3); hold on;
semilogy(Dl, max(Md, 1e-17), 'o', 'MarkerFaceColor', 'auto');
set(gca, 'XScale', 'log'); xticks(Dl); grid on; ylim([1e-17 1]);
xlabel('D'); ylabel('\epsilon'); title('Gauss x sinc \alpha=4, half-span 4 T_s'); legend('prediction', 'measured', 'Location', 'east');
subplot(1,3,3);
allP = [Eprd(:); Pd(:)];  allM = [Emea(:); Md(:)];
loglog(allP, allM, 'o'); hold on; plot([1e-4 1], [1e-4 1], 'k--');
axis([1e-4 1 1e-4 1]); grid on;
xlabel('predicted \epsilon (from R[kD])'); ylabel('measured \epsilon'); title('prediction vs measurement');
ok = allP > 1e-6;
fprintf('max relative deviation, measured vs predicted: %.1f %%\n', 100*max(abs(allM(ok) - allP(ok))./allP(ok)));

%% ---- Part F: what interference does to a constellation (Figure 6) -------
M2 = 32;  N2 = 16;  L2 = M2*N2*D;
X2 = ((2*randi(4,M2,N2)-5) + 1j*(2*randi(4,M2,N2)-5))/sqrt(10);
[g1, e1] = pulse_rect(D);
[g2, e2] = pulse_gs(D, 33, 4, M2);
[g3, e3] = pulse_orth(D, 65, 16, M2, L2);
PG = {g1, g2, g3};  PE = [e1 e2 e3];
pn = {'rect', 'Gauss x sinc, \alpha=4, 33 taps', 'orthogonalised, \alpha=16, 65 taps'};
figure('Name','Fig6: constellations after loopback');
for i = 1:3
    Y2 = zak_sampled(izak_sampled(X2, PG{i}, PE(i), D), M2, N2, PG{i}, PE(i), D);
    evm = 100*norm(Y2(:) - X2(:))/norm(X2(:));
    [~,~,~,~,ep] = pulse_stats(PG{i}, D);
    fprintf('%-34s EVM = %.4f %%   (100*eps_pred = %.4f %%)\n', pn{i}, evm, 100*ep);
    subplot(1,3,i);
    plot(real(Y2(:)), imag(Y2(:)), '.', 'MarkerSize', 6); hold on;
    plot(real(X2(:)), imag(X2(:)), 'r+', 'MarkerSize', 8);
    axis equal; grid on; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
    xlabel('Re'); ylabel('Im'); title(sprintf('%s: EVM %.3g %%', pn{i}, evm));
end
sgtitle('blue = Y after IZak + Zak (no channel, no noise), red + = transmitted X');

%% ---- Part G: why the pulse is normalised (slides p3, p6) (Figure 7) -----
Wv2 = 1:8;  cvals = zeros(size(Wv2));  gain = zeros(size(Wv2));
for j = 1:numel(Wv2)
    Lgj = 2*D*Wv2(j) + 1;
    ellj = (-(Lgj-1)/2 : (Lgj-1)/2).';
    graw = exp(-pi*4*((ellj/D)/M).^2) .* sincn(ellj/D);          % NOT normalised
    cvals(j) = sum(graw.^2);
    sraw = izak_sampled(X, graw, -(Lgj-1)/2, D);
    Yraw = zak_sampled(sraw, M, N, graw, -(Lgj-1)/2, D);
    gain(j) = sum(real(conj(X(:)).*Yraw(:)))/sum(abs(X(:)).^2);   % least-squares gain Y ~ gain*X
end
disp('half-span | sum g^2 | measured gain | difference (%)');
disp([Wv2.' cvals.' gain.' 100*abs(gain.'-cvals.')./cvals.']);
figure('Name','Fig7: normalisation');
subplot(1,2,1);
plot(Wv2, cvals, 'o-', 'LineWidth', 1.5); grid on;
xlabel('pulse half-span (T_s)'); ylabel('c = \Sigma g^2 (not normalised)'); title('energy of the raw pulse');
subplot(1,2,2);
plot(cvals, gain, 'o'); hold on; plot([3 4], [3 4], 'k--'); grid on; axis equal;
xlabel('c = \Sigma g^2'); ylabel('measured gain of the loopback'); title('without normalisation Y \approx c X');

%% ==================== local functions (must stay at the end) ==============
function [g, ell0] = pulse_rect(D)
    g = ones(D,1)/sqrt(D);  ell0 = 0;                  % unit-energy rectangle, width Ts
end

function y = sincn(x)
    y = ones(size(x));  nz = x ~= 0;
    y(nz) = sin(pi*x(nz))./(pi*x(nz));                 % sin(pi x)/(pi x), 1 at x = 0
end

function [g, ell0] = pulse_gs(D, Lg, alpha, M)
% Gaussian x sinc (slide p15), sampled at Tsamp = Ts/D, Lg (odd) taps centred on 0,
% tau_s = Ts, T = M*Ts, normalised to unit energy (slide p3).
    c   = (Lg-1)/2;
    ell = (-c:c).';
    t   = ell/D;                                       % time in units of Ts
    g   = exp(-pi*alpha*(t/M).^2) .* sincn(t);
    g   = g / sqrt(sum(abs(g).^2));
    ell0 = -c;
end

function [g, ell0] = pulse_orth(D, Lg, alpha, M, L)
% Orthogonalised Gauss x sinc: make the power spectrum folded by D flat, so that
% sum_a g[a] g[a - kD] = delta_k exactly (circular frame of length L), then truncate to Lg taps.
    Lb = 32*D + 1;  cb = (Lb-1)/2;
    ell  = (-cb:cb).';  t = ell/D;
    base = exp(-pi*alpha*(t/M).^2) .* sincn(t);
    base = base / sqrt(sum(base.^2));
    full = zeros(L, 1);
    full(mod(ell, L) + 1) = base;                      % pulse centred at index 0 (circular)
    Gf = fft(full);
    Lq = L/D;
    F  = zeros(Lq, 1);
    for k = 0:D-1
        F = F + abs(Gf(k*Lq + (1:Lq))).^2;            % fold the power spectrum by D
    end
    Go = Gf .* sqrt(D ./ repmat(F, D, 1));
    go = real(ifft(Go));
    cg = (Lg-1)/2;
    g  = go(mod((-cg:cg).', L) + 1);
    g  = g / sqrt(sum(g.^2));
    ell0 = -cg;
end

function [lags, Rfull, kk, Rk, ep] = pulse_stats(g, D)
% Autocorrelation of the pulse, its samples at multiples of D (the Gram taps G_k),
% and eps = sqrt(sum_k |G_k - delta_k|^2).
    Lg    = numel(g);
    Rfull = conv(g, flipud(conj(g)));                  % Rfull(i) <-> lag = i - Lg
    lags  = (-(Lg-1):(Lg-1)).';
    sel   = mod(lags, D) == 0;
    kk    = lags(sel)/D;
    Rk    = Rfull(sel);
    ep    = sqrt(sum(abs(Rk - (kk == 0)).^2));
end

function H = build_H(g, ell0, D, M, N)
% End-to-end operator: column i = loopback response to a unit symbol at DD cell i (i = m+1 + M n).
    H = zeros(M*N, M*N);
    for i = 1:M*N
        Xi = zeros(M, N);  Xi(i) = 1;
        s  = izak_sampled(Xi, g, ell0, D);
        Y  = zak_sampled(s, M, N, g, ell0, D);
        H(:, i) = Y(:);
    end
end

function e = meas_eps(g, ell0, D, M, N, T)
% rms relative loopback error over T random 16-QAM-like frames
    acc = 0;
    for t = 1:T
        Xr = ((2*randi(4,M,N)-5) + 1j*(2*randi(4,M,N)-5))/sqrt(10);
        Y  = zak_sampled(izak_sampled(Xr, g, ell0, D), M, N, g, ell0, D);
        acc = acc + sum(abs(Y(:) - Xr(:)).^2)/sum(abs(Xr(:)).^2);
    end
    e = sqrt(acc/T);
end

function s = izak_sampled(X, g, ell0, D)
% s[a] = sum_m sum_q v_m[q] g[a - D m - P q],  v_m[q] = (1/N) sum_n X[m,n] e^{+j2pi nq/N}
    [M, N] = size(X);
    P = M*D;  L = N*P;
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
% u_m[q] = sum_ell r[D m + P q + ell] conj(g[ell]);  Y[m,n] = sum_q u_m[q] e^{-j2pi nq/N}
    L = M*N*D;  r = r(:);
    acc = zeros(L, 1);
    for i = 1:numel(g)
        acc = acc + circshift(r, -(ell0 + i - 1)) * conj(g(i));
    end
    k = (0:M-1).' + M*(0:N-1);
    U = reshape(acc(D*k(:) + 1), M, N);
    Y = fft(U, [], 2);
end