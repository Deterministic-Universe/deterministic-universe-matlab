% Requires MATLAB R2016b or newer (local functions at the end of a script).

% --- 1. The Physics Engine: 2D incompressible Navier-Stokes in vorticity form ---
% Taking the curl of the momentum equation removes the pressure:
%     d(w)/dt + u . grad(w) = nu * laplacian(w),      w = dv/dx - du/dy,      nu = mu / rho
% Incompressibility (div u = 0) is built in through a streamfunction:  u = dpsi/dy, v = -dpsi/dx, laplacian(psi) = -w.
N = 128; nu = 5e-4; dt = 0.002; t_end = 1.2;
x = (0:N-1) / N;
[X, Y] = meshgrid(x, x);                           % periodic unit square
k = 2 * pi * [0:N/2-1, -N/2:-1];                   % spectral wavenumbers (FFT ordering)
[KX, KY] = meshgrid(k, k);
K2 = KX.^2 + KY.^2;
K2s = K2; K2s(1, 1) = 1;                           % avoid dividing by zero at the mean mode
mask = (abs(KX) < (2/3) * max(k)) & (abs(KY) < (2/3) * max(k));   % 2/3 rule against aliasing
cfg = struct('KX', KX, 'KY', KY, 'K2', K2, 'K2s', K2s, 'mask', mask, 'nu', nu);

% --- 2. Initial Conditions: a double shear layer with a small perturbation (it rolls up into vortices) ---
rho_s = 30; delta = 0.05;
u0 = (Y <= 0.5) .* tanh(rho_s * (Y - 0.25)) + (Y > 0.5) .* tanh(rho_s * (0.75 - Y));
v0 = delta * sin(2 * pi * X);
w_hat = (1i * KX .* fft2(v0) - 1i * KY .* fft2(u0)) .* mask;
w_start = real(ifft2(w_hat));

% --- 3. Time stepping: classical Runge-Kutta 4 ---
[u, v] = velocity(w_hat, cfg);
E0 = 0.5 * mean(u(:).^2 + v(:).^2);
for step = 1:round(t_end / dt)
    k1 = dwdt(w_hat, cfg);
    k2 = dwdt(w_hat + 0.5 * dt * k1, cfg);
    k3 = dwdt(w_hat + 0.5 * dt * k2, cfg);
    k4 = dwdt(w_hat + dt * k3, cfg);
    w_hat = w_hat + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
[u, v] = velocity(w_hat, cfg);
fprintf('kinetic energy fraction remaining: %.4f (viscosity removes energy)\n', 0.5 * mean(u(:).^2 + v(:).^2) / E0);
div = real(ifft2(1i * KX .* fft2(u) + 1i * KY .* fft2(v)));
fprintf('max |div u| (should be ~0): %.3e\n', max(abs(div(:))));

% --- 4. Generative Visual Mapping: vorticity ---
w_end = real(ifft2(w_hat));
lim = max(abs(w_end(:)));
nc = 128;                                          % a blue-white-red colormap for signed vorticity
diverging = [linspace(0, 1, nc)', linspace(0, 1, nc)', ones(nc, 1); ...
             ones(nc, 1), linspace(1, 0, nc)', linspace(1, 0, nc)'];
figure('Color', 'k', 'Position', [100 100 1100 550]);
fields = {w_start, w_end};
titles = {'Vorticity, t = 0', sprintf('Vorticity, t = %g', t_end)};
for i = 1:2
    subplot(1, 2, i);
    imagesc(x, x, fields{i}); set(gca, 'YDir', 'normal');
    caxis([-lim lim]); colormap(gca, diverging); axis equal off;
    title(titles{i}, 'Color', 'w');
end

% --- Local functions (they go at the end of a script) ---
function [u, v] = velocity(w_hat, cfg)
    psi_hat = w_hat ./ cfg.K2s;
    psi_hat(1, 1) = 0;
    u = real(ifft2( 1i * cfg.KY .* psi_hat));      % u =  dpsi/dy
    v = real(ifft2(-1i * cfg.KX .* psi_hat));      % v = -dpsi/dx
end

function d = dwdt(w_hat, cfg)
    [u, v] = velocity(w_hat, cfg);
    wx = real(ifft2(1i * cfg.KX .* w_hat));
    wy = real(ifft2(1i * cfg.KY .* w_hat));
    advection = fft2(u .* wx + v .* wy) .* cfg.mask;   % nonlinear term  u . grad(w)
    d = -advection - cfg.nu * cfg.K2 .* w_hat;         % advection + viscous diffusion
end
