%% Chapter 14 -- Electromagnetism and Classical Field Theory
% The Deterministic Universe -- Anthony Briggs
% Requires MATLAB R2016b or newer. Local functions are at the end of the file.
%
% Part A  Lorentz force on a charge:  m a = q (E + v x B)
%         uniform E and B give a cycloid with E x B drift  v_d = E x B / B^2
% Part B  Maxwell's equations on a 2-D Yee grid (TM mode: Ez, Hx, Hy), eps0 = mu0 = c = 1
%         dHx/dt = -dEz/dy      dHy/dt = +dEz/dx      dEz/dt = dHy/dx - dHx/dy - Jz
%         field energy  W = sum( 1/2 Ez^2 + 1/2 Hx^2 + 1/2 Hy^2 ) dx^2
%         Poynting theorem: dW/dt = - sum(J.E) dx^2 - (flux through the boundary)
clear; clc; close all;

%% 1 Physics Engine -- Part A: charged particle
q = 1; m = 1; Ey = 0.2; Bz = 1.0;             % E = (0, Ey), B = (0, 0, Bz)
lorentz = @(t, s) [ s(3);
                    s(4);
                    (q/m) * ( s(4) * Bz);
                    (q/m) * ( Ey - s(3) * Bz) ];

%% 2 Initial Conditions -- Part A
sA0   = [0; 0; 0; 1.0];                        % x, y, vx, vy
tA    = linspace(0, 8*pi, 1601);
opts  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, SA] = ode45(lorentz, tA, sA0, opts);
x = SA(:,1); y = SA(:,2); vx = SA(:,3); vy = SA(:,4);
EA = 0.5*(vx.^2 + vy.^2) - Ey*y;                % K + U with U = -q E y

%% 1/2 Part B: Yee-grid FDTD
N = 201; dx = 1.0; dt = 0.5;                   % Courant number c dt/dx = 0.5 (< 1/sqrt(2))
i0 = (N + 1)/2; j0 = i0;                       % source cell (centre)
Ez = zeros(N, N); Hx = zeros(N, N-1); Hy = zeros(N-1, N);
om = 0.5; Tb = 50;                             % drive frequency and burst length
Jfun = @(t) (t >= 0 && t <= Tb) * sin(om*t) * sin(pi*t/Tb)^2;
nsteps = round(70/dt);
tt = zeros(1, nsteps); Wf = zeros(1, nsteps); Wk = zeros(1, nsteps);
work = 0; snap = [];
for k = 1:nsteps
    t = (k-1)*dt;
    Hxo = Hx; Hyo = Hy;
    Hx = Hx - dt/dx * (Ez(:,2:end) - Ez(:,1:end-1));
    Hy = Hy + dt/dx * (Ez(2:end,:) - Ez(1:end-1,:));
    Jn = Jfun(t + dt/2);
    Ezold = Ez(i0, j0);
    Ez(2:N-1, 2:N-1) = Ez(2:N-1, 2:N-1) + dt/dx * ( ...
        (Hy(2:end, 2:N-1) - Hy(1:end-1, 2:N-1)) - ...
        (Hx(2:N-1, 2:end) - Hx(2:N-1, 1:end-1)) );
    Ez(i0, j0) = Ez(i0, j0) - dt*Jn;
    work = work - Jn * 0.5*(Ezold + Ez(i0, j0)) * dx*dx * dt;   % energy the source gives the field
    Wf(k) = 0.5*sum(Ez(:).^2) + 0.5*sum(Hxo(:).*Hx(:)) + 0.5*sum(Hyo(:).*Hy(:));
    Wk(k) = work; tt(k) = k*dt;
    if k == round(60/dt), snap = Ez; end
end

%% 3 Generative Visual Mapping + ledger
fprintf('--- Checking your run (Chapter 14) ---\n');
fprintf('Part A: E = K - q E y: start %.4f, max drift %.2e\n', EA(1), max(abs(EA - EA(1))));
fprintf('Part A: x after one gyro-period = %.4f  (drift*2pi = %.4f)\n', interp1(tA, x, 2*pi), 0.2*2*pi);
fprintf('Part A: y range [%.4f, %.4f]\n', min(y), max(y));
for tq = [40 50 60 70]
    kk = round(tq/dt);
    fprintf('Part B: t=%2d  W = %.5f   work by source = %.5f\n', tq, Wf(kk), Wk(kk));
end
kk = round(55/dt);
fprintf('Part B: spread of W after the source is off = %.2e\n', max(Wf(kk:end)) - min(Wf(kk:end)));
row = snap(i0+15:i0+49, j0);
zc  = find(row(1:end-1).*row(2:end) < 0);
fprintf('Part B: wavelength from zero crossings = %.1f  (theory 2*pi*c/om = %.2f)\n', 2*mean(diff(zc)), 2*pi/om);

fig = figure('Color', 'k', 'Position', [100 100 1500 460]);
ax1 = subplot(1, 3, 1);
plot(x, y, 'Color', [0 0.9 1], 'LineWidth', 1.6); axis equal;
title('Lorentz force: E x B drift', 'Color', 'w'); xlabel('x'); ylabel('y'); dark_axes(ax1);

ax2 = subplot(1, 3, 2);
mm = 0.6*max(abs(snap(:)));
imagesc(snap.', [-mm mm]); axis xy; colormap(ax2, redblue(256));
title('E_z at t = 60', 'Color', 'w'); xlabel('i'); ylabel('j'); dark_axes(ax2);

ax3 = subplot(1, 3, 3); hold on;
plot(tt, Wf, 'Color', [0 1 0.53], 'LineWidth', 1.5);
plot(tt, Wk, '--', 'Color', [1 0.9 0], 'LineWidth', 1.5);
title('Ledger', 'Color', 'w'); xlabel('t');
lg = legend('field energy W', 'work done by source', 'Location', 'southeast'); style_legend(lg); dark_axes(ax3);

set(fig, 'InvertHardcopy', 'off');
print(fig, 'fig16_14.png', '-dpng', '-r300');

%% ---------------- local functions ----------------
function dark_axes(ax)
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'Box', 'on');
end

function style_legend(lg)
    set(lg, 'Color', 'k', 'EdgeColor', 'w');
    try, set(lg, 'TextColor', 'w'); catch, end   % TextColor needs R2018a+
end

function cmap = redblue(n)
    % blue -> white -> red diverging map, no toolbox needed
    h = floor(n/2);
    up   = linspace(0, 1, h)';
    down = linspace(1, 0, n - h)';
    cmap = [ [up, up, ones(h,1)]; [ones(n-h,1), down, down] ];
end
