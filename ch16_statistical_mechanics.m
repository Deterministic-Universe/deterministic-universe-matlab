%% Chapter 16 -- Thermodynamics and Statistical Mechanics: a Hard-Disk Gas
% The Deterministic Universe -- Anthony Briggs
% Requires MATLAB R2016b or newer (implicit expansion). Local functions at the end.
%
% Dynamics   free flight, elastic wall bounces, elastic disk-disk collisions
%            (equal masses: the velocity components along the line of centres swap)
% Boltzmann  S = k ln(Omega); here the coarse-grained form S = -sum p_i ln p_i
%            over an M x M grid of cells (k_B = 1)
% Maxwell-Boltzmann (2-D)  f(v) = (m v / kT) exp(-m v^2 / 2kT)
%            mean speed sqrt(pi kT/2), rms speed sqrt(2 kT), most probable speed sqrt(kT)
clear; clc; close all;
rng(7);

%% 1 Physics Engine
N = 1000; r = 0.005; L = 1.0; kT = 1.0; m = 1;
dt = 0.002; nsteps = 1500; M = 10;

%% 2 Initial Conditions -- every disk starts in the LEFT half of the box
pos = zeros(N, 2); n = 0;
while n < N
    p = [r + (0.5 - 2*r)*rand, r + (L - 2*r)*rand];
    if n == 0 || min(hypot(pos(1:n,1) - p(1), pos(1:n,2) - p(2))) > 2*r
        n = n + 1; pos(n,:) = p;
    end
end
v = randn(N, 2);
v = v - mean(v, 1);                           % zero total momentum
v = v * sqrt(2*N*kT / sum(v(:).^2));          % scale so that K = N kT exactly
K0 = 0.5 * sum(v(:).^2);

t = zeros(1, nsteps+1); S = t; K = t;
S(1) = cell_entropy(pos, L, M); K(1) = K0;
ncoll = 0; nlate = 0;
for s = 1:nsteps
    pos = pos + v*dt;
    for dcol = 1:2                            % elastic reflection from the four walls
        lo = pos(:,dcol) < r;      hi = pos(:,dcol) > L - r;
        pos(lo,dcol) = 2*r - pos(lo,dcol);        v(lo,dcol) = -v(lo,dcol);
        pos(hi,dcol) = 2*(L - r) - pos(hi,dcol);  v(hi,dcol) = -v(hi,dcol);
    end
    dxm = pos(:,1) - pos(:,1)';  dym = pos(:,2) - pos(:,2)';
    [I, J] = find(triu(dxm.^2 + dym.^2 < (2*r)^2, 1));
    for q = 1:numel(I)
        i = I(q); j = J(q);
        dr = pos(i,:) - pos(j,:);  dist = norm(dr);
        if dist == 0, continue; end
        nv = dr / dist;  vn = (v(i,:) - v(j,:)) * nv';
        if vn < 0                              % approaching: exchange normal components
            v(i,:) = v(i,:) - vn*nv;  v(j,:) = v(j,:) + vn*nv;
            ncoll = ncoll + 1;
            if s*dt > 1.0, nlate = nlate + 1; end
        end
    end
    t(s+1) = s*dt; S(s+1) = cell_entropy(pos, L, M); K(s+1) = 0.5*sum(v(:).^2);
end

%% 3 Generative Visual Mapping + ledger
sp = hypot(v(:,1), v(:,2));
fprintf('--- Checking your run (Chapter 16) ---\n');
fprintf('K0 = %.1f   max |K-K0| = %.2e\n', K0, max(abs(K - K0)));
fprintf('S(0) = %.3f  (theory ln(M^2/2) = %.3f)\n', S(1), log(M^2/2));
fprintf('S(end) = %.3f  (ceiling ln(M^2) = %.3f ; finite-N offset (M^2-1)/(2N) = %.3f)\n', ...
        S(end), log(M^2), (M^2 - 1)/(2*N));
for tq = [0.25 0.5 1 2 3]
    fprintf('t = %.2f   S = %.3f\n', tq, S(round(tq/dt) + 1));
end
fprintf('mean speed = %.3f (theory %.3f)   rms speed = %.3f (theory %.3f)\n', ...
        mean(sp), sqrt(pi*kT/2), sqrt(mean(sp.^2)), sqrt(2*kT));
fprintf('collisions per disk per unit time (t>1) = %.1f   (dilute theory %.1f)\n', ...
        2*nlate/N/(t(end) - 1), 2*(2*r)*(N/L^2)*sqrt(2)*sqrt(pi*kT/2));

col = {[0 0.9 1], [1 0 1], [1 0.9 0], [0 1 0.53], [1 0.53 0]};
fig = figure('Color', 'k', 'Position', [100 100 1500 460]);
ax1 = subplot(1, 3, 1);
plot(pos(:,1), pos(:,2), '.', 'Color', col{1}, 'MarkerSize', 6); axis equal; xlim([0 1]); ylim([0 1]);
title('1000 disks at t = 3.0', 'Color', 'w'); dark_axes(ax1);

ax2 = subplot(1, 3, 2); hold on;
histogram(sp, 30, 'Normalization', 'pdf', 'FaceColor', col{2}, 'EdgeColor', 'k');
vv = linspace(0, 4, 200);
plot(vv, (m*vv/kT) .* exp(-m*vv.^2/(2*kT)), 'Color', col{3}, 'LineWidth', 1.8);
title('Speed distribution', 'Color', 'w'); xlabel('speed');
lg = legend('simulation', 'Maxwell-Boltzmann (2D)'); style_legend(lg); dark_axes(ax2);

ax3 = subplot(1, 3, 3); hold on;
plot(t, S, 'Color', col{4}, 'LineWidth', 1.6);
plot([0 t(end)], log(M^2)*[1 1], '--', 'Color', [0.6 0.6 0.6]);
plot(t, K/K0, 'Color', col{5}, 'LineWidth', 1.6);
title('Ledger', 'Color', 'w'); xlabel('t');
lg = legend('S (coarse-grained)', 'ln(M^2)', 'K / K_0', 'Location', 'east'); style_legend(lg); dark_axes(ax3);

set(fig, 'InvertHardcopy', 'off');
print(fig, 'fig16_16.png', '-dpng', '-r300');

%% ---------------- local functions ----------------
function S = cell_entropy(pos, L, M)
    ix = min(floor(pos(:,1)/L*M), M-1);
    iy = min(floor(pos(:,2)/L*M), M-1);
    cnt = accumarray(ix*M + iy + 1, 1, [M*M, 1]);
    p = cnt(cnt > 0) / size(pos, 1);
    S = -sum(p .* log(p));
end

function dark_axes(ax)
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'Box', 'on');
end

function style_legend(lg)
    set(lg, 'Color', 'k', 'EdgeColor', 'w');
    try, set(lg, 'TextColor', 'w'); catch, end   % TextColor needs R2018a+
end
