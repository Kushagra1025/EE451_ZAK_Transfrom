clear; clc; close all;
M = 128;  N = 128;  B = 500e6;

%Part A....

for D = [1 4]
    p = otfs_params(M, N, B, D);
    fprintf('---- M=%d, N=%d, B=%.0f MHz, D=%d ----\n', M, N, B/1e6, D);
    fprintf('Ts (delay bin)        = %.3f ns\n',  p.Ts*1e9);
    fprintf('T = M*Ts (Zak period) = %.3f ns\n',  p.T*1e9);
    fprintf('Tframe = N*T          = %.4f us\n',  p.Tframe*1e6);
    fprintf('Tsamp = Ts/D          = %.3f ns (fs = %.0f MS/s)\n', p.Tsamp*1e9, p.fs/1e6);
    fprintf('P = M*D               = %d samples per Zak period\n', p.P);
    fprintf('Samples per frame     = %d\n', p.Lframe);
    fprintf('Frame rate            = %.3f frames/s\n', p.frameRate);
    fprintf('Doppler spacing       = %.4f kHz\n', p.dnu/1e3);
    fprintf('Doppler span (1/T)    = %.5f MHz\n\n', p.nuSpan/1e6);
    assert(abs(p.frameRate - B/(M*N)) < 1e-6);   % slide p19: f_frame = B/(MN)
    assert(p.Lframe == M*N*D);
end

%Part B....

M = 8;  N = 4;  D = 2;  P = M*D;
[mm, qq] = ndgrid(0:M-1, 0:N-1);     % M x N grids of m and q
A = D*mm + P*qq;                     % sample index of each pulse centre

assert(numel(unique(A)) == M*N);                        % all positions distinct
assert(isequal(sort(A(:)).', D*(0:M*N-1)));             % they tile the timeline in steps of D

figure('Name','Index map');
imagesc(0:N-1, 0:M-1, A); axis xy; colormap(summer); colorbar;
xlabel('Zak period index q'); ylabel('Delay bin m');
title(sprintf('Sample index a = D m + P q   (M=%d, N=%d, D=%d)', M, N, D));
for m = 0:M-1
    for q = 0:N-1
        text(q, m, num2str(A(m+1,q+1)), 'HorizontalAlignment','center');
    end
end
%saveas(gcf, 'S1_fig1_index_map.png');

%Part C....

M = 8;  N = 4;  D = 4;  P = M*D;  Lframe = N*P;     % 128 samples

figure('Name','Timeline'); hold on;
plot(0:Lframe-1, zeros(1,Lframe), 'k.', 'MarkerSize', 6);   % every ADC/DAC sample
cols = lines(N);
for q = 0:N-1
    c = D*(0:M-1) + P*q;                                    % pulse centres in period q
    stem(c, ones(1,M), 'Color', cols(q+1,:), 'MarkerFaceColor', cols(q+1,:), 'LineWidth', 1.2);
    text(P*q + P/2, 1.25, sprintf('period q=%d', q), 'HorizontalAlignment','center');
end
for q = 0:N, xline(P*q, 'k--'); end                         % Zak-period boundaries
xlim([-2 Lframe+2]); ylim([-0.2 1.5]); yticks([]);
xlabel('ADC/DAC sample index a');
title(sprintf('Pulse centres (M=%d, N=%d, D=%d, P=%d): one frame = %d samples', M, N, D, P, Lframe));
%saveas(gcf, 'S1_fig2_timeline.png');

%Part D....

M = 16;  N = 8;  rng(1);
X = zeros(M, N);
qpsk = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
K = 10;
idx = randperm(M*N, K);
X(idx) = qpsk(randi(4, 1, K)) .* (0.5 + rand(1, K));   % random amplitudes, only for visibility
[mm, nn] = ndgrid(0:M-1, 0:N-1);

figure('Name','DD grid');
subplot(1,2,1);
imagesc(0:N-1, 0:M-1, abs(X)); axis xy; colormap(gca, flipud(gray)); colorbar;
xlabel('Doppler bin n'); ylabel('Delay bin m'); title('|X[m,n]| heat map');
subplot(1,2,2);
stem3(nn(idx), mm(idx), abs(X(idx)), 'filled'); grid on;
xlim([0 N-1]); ylim([0 M-1]);
xlabel('Doppler bin n'); ylabel('Delay bin m'); zlabel('|X|'); title('|X[m,n]| stems');
%saveas(gcf, 'S1_fig3_DD_grid.png');

%Part E....

figure('Name','Resolutions');
tiledlayout(1,2);
nexttile;
Bv = logspace(7, 9.5, 100);
loglog(Bv/1e6, 1e9./Bv, 'LineWidth', 1.5); grid on;
xlabel('Bandwidth B (MHz)'); ylabel('\Delta\tau = T_s = 1/B (ns)');
title('Delay resolution');
nexttile;
Nv = 2.^(3:10);  Mfix = 128;  Bfix = 500e6;
dnu = 1 ./ (Nv*Mfix/Bfix);
loglog(Nv, dnu/1e3, 'o-', 'LineWidth', 1.5); grid on;
xlabel('N (periods per frame)'); ylabel('\Delta\nu = 1/T_{frame} (kHz)');
title('Doppler spacing (M=128, B=500 MHz)');
%saveas(gcf, 'S1_fig4_resolutions.png');