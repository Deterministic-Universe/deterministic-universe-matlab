% --- 1. Resonance parameters ---
m_c = 5; n_c = 4;                      % circle: m nodal diameters, n-th radial mode (n - 1 interior nodal circles)
m_s = 2; n_s = 5;                      % square plate mode numbers
res = 700;

% --- 2. The Physics Engines (modal solutions) ---
% Circular membrane (a drumhead) of radius 1:  u = J_m(alpha_mn r) cos(m theta)
% MATLAB has no built-in table of Bessel zeros, so find alpha_mn = n-th positive zero of J_m
% by scanning for sign changes of J_m and polishing each with fzero.
xg = linspace(0.1, 100, 40000);        % the scan range is enough for m, n up to about 20
sg = besselj(m_c, xg);
idx = find(sg(1:end-1) .* sg(2:end) < 0);
alpha = fzero(@(x) besselj(m_c, x), [xg(idx(n_c)), xg(idx(n_c) + 1)]);
fprintf('alpha(%d,%d) = %.6f\n', m_c, n_c, alpha);

x = linspace(-1, 1, res);
[X1, Y1] = meshgrid(x, x);
R = hypot(X1, Y1);
Theta = atan2(Y1, X1);
U1 = besselj(m_c, alpha * R) .* cos(m_c * Theta);
U1(R > 1) = NaN;                       % outside the drum

% Idealized square plate:  u = cos(n pi x) cos(m pi y) - cos(m pi x) cos(n pi y)
xs = linspace(0, 1, res);
[X2, Y2] = meshgrid(xs, xs);
U2 = cos(n_s * pi * X2) .* cos(m_s * pi * Y2) - cos(m_s * pi * X2) .* cos(n_s * pi * Y2);

% --- 3. Generative Visual Mapping ---
nc = 128;                              % a blue-white-red colormap for signed amplitude
diverging = [linspace(0, 1, nc)', linspace(0, 1, nc)', ones(nc, 1); ...
             ones(nc, 1), linspace(1, 0, nc)', linspace(1, 0, nc)'];
sandmap = [linspace(1, 0.5, 64)', linspace(1, 0.2, 64)', linspace(0.9, 0, 64)'];

figure('Color', 'k', 'Position', [100 100 1500 500]);

% 'Mandala': amplitude with the nodal lines (u = 0) drawn in white
subplot(1, 3, 1);
peak = max(abs(U1(:)));
contourf(X1, Y1, U1, linspace(-peak, peak, 101), 'LineStyle', 'none');
hold on; contour(X1, Y1, U1, [0 0], 'w', 'LineWidth', 0.7);
colormap(gca, diverging); axis equal off;
title(sprintf('Mandala: drumhead mode (m=%d, n=%d)', m_c, n_c), 'Color', 'w');

% 'Sand': sand collects where the surface stays still (u is close to 0)
subplot(1, 3, 2);
sand = exp(-(U2 / (0.06 * max(abs(U2(:))))).^2);
imagesc(xs, xs, sand); set(gca, 'YDir', 'normal');
colormap(gca, sandmap); axis equal off;
title(sprintf('Sand figure: square plate (m=%d, n=%d)', m_s, n_s), 'Color', 'w');

subplot(1, 3, 3);
peak2 = max(abs(U2(:)));
contourf(X2, Y2, U2, linspace(-peak2, peak2, 101), 'LineStyle', 'none');
hold on; contour(X2, Y2, U2, [0 0], 'w', 'LineWidth', 0.7);
colormap(gca, diverging); axis equal off;
title('Same plate: amplitude and nodal lines', 'Color', 'w');
