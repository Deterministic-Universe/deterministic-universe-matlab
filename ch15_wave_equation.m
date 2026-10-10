%% Chapter 15 -- The Wave Equation: FDTD double-slit diffraction
% Self-contained: paste into the editor or the command window and run.
% No local functions, no other files needed. Needs MATLAB R2018a or newer.
%   u_tt = c^2 (u_xx + u_yy)
%   leapfrog: u(n+1) = 2u(n) - u(n-1) + (c dt/dx)^2 * Laplacian(u(n)),  c dt/dx <= 1/sqrt(2)
%   energy:   E = sum(1/2 u_t^2 + 1/2 c^2 |grad u|^2)
clear; clc; close all;

%% 1 Physics Engine
Nx = 360; Ny = 400; dx = 1.0; c = 1.0; dt = 0.5;
lam = 10; kw = 2*pi/lam;                         % carrier wavelength, wavenumber
x0 = 60; sg = 15;                                % pulse centre and envelope width
xw = 150; a = 6; d = 40; cy = Ny/2; xs = 300;    % wall, slit width, slit spacing, centre line, screen
r2 = (c*dt/dx)^2;
[X, Y] = ndgrid((0:Nx-1)*dx, (0:Ny-1)*dx);       % X = x, Y = y (0-based coordinates)
pulse = @(xx) exp(-(xx - x0).^2/(2*sg^2)) .* cos(kw*(xx - x0));
gd = @(u1,u2) sum(sum((u1(2:end,:)-u1(1:end-1,:)).*(u2(2:end,:)-u2(1:end-1,:)))) ...
            + sum(sum((u1(:,2:end)-u1(:,1:end-1)).*(u2(:,2:end)-u2(:,1:end-1))));

%% 2 Initial Conditions
wall = abs(X - xw) <= 1;                         % wall three cells thick
for yc = [cy - d/2, cy + d/2]
    wall(abs(X - xw) <= 1 & Y >= yc - a/2 & Y < yc + a/2) = false;   % open the slits
end
mask = ~wall;
mask(1,:) = false; mask(end,:) = false; mask(:,1) = false; mask(:,end) = false;
un = pulse(X) .* mask;                           % u at step n
up = pulse(X + c*dt) .* mask;                    % u at step n-1

nsteps = round(330/dt);
tt = zeros(1,nsteps); Ek = tt; Es = tt;
Iscreen = zeros(Ny,1); snap = zeros(Nx,Ny);
for k = 1:nsteps
    lap = zeros(Nx,Ny);
    lap(2:end-1,2:end-1) = un(3:end,2:end-1) + un(1:end-2,2:end-1) ...
                         + un(2:end-1,3:end) + un(2:end-1,1:end-2) - 4*un(2:end-1,2:end-1);
    unew = (2*un - up + r2*lap) .* mask;
    Ek(k) = 0.5*sum(sum(((unew - un)/dt).^2))*dx*dx;
    Es(k) = 0.5*c^2*gd(unew, un);
    tt(k) = k*dt;
    up = un; un = unew;
    Iscreen = Iscreen + (un(xs+1,:).^2).'*dt;
    if k == round(200/dt), snap = un; end
end
Et = Ek + Es;

%% 3 Checks and figure
fprintf('--- Checking your run (Chapter 15) ---\n');
fprintf('E total = %.2f, relative drift = %.2e\n', Et(1), max(abs(Et-Et(1)))/Et(1));
for tq = [10 100 150]
    kk = round(tq/dt);
    fprintf('t = %3d: kinetic = %7.1f  strain = %7.1f\n', tq, Ek(kk), Es(kk));
end
Inorm = Iscreen/max(Iscreen);
yy = (0:Ny-1)' - cy;
isp = find(Inorm(2:end-1) > Inorm(1:end-2) & Inorm(2:end-1) >= Inorm(3:end) & Inorm(2:end-1) > 0.3) + 1;
fprintf('fringe maxima at y - yc = %s\n', mat2str(yy(isp)'));
fprintf('theory: spacing lam*L/d = %.1f; exact first order = %.1f\n', ...
        lam*(xs-xw)/d, (xs-xw)*tan(asin(lam/d)));

dk  = @(ax) set(ax,'Color','k','XColor','w','YColor','w','Box','on');
lgs = @(lg) set(lg,'Color','k','EdgeColor','w','TextColor','w');
h = 128;
rb = [linspace(0,1,h)' linspace(0,1,h)' ones(h,1); ones(h,1) linspace(1,0,h)' linspace(1,0,h)'];
fig = figure('Color','k','Position',[100 100 1500 480]);

ax = subplot(1,3,1);
mm = max(abs(snap(:)));
imagesc(0:Nx-1, 0:Ny-1, snap.', [-mm mm]); axis xy; colormap(ax, rb);
hold on; plot([xw xw], [0 Ny], 'w-', 'LineWidth', 0.8);
title('u at t = 200 (double slit)','Color','w'); xlabel('x'); ylabel('y'); dk(ax);

ax = subplot(1,3,2);
plot(yy, Inorm, 'Color', [0 0.9 1], 'LineWidth', 1.5); xlim([-130 130]);
title('Intensity on the screen (x = 300)','Color','w'); xlabel('y - y_c');
ylabel('\int u^2 dt (normalised)'); dk(ax);

ax = subplot(1,3,3); hold on;
plot(tt, Ek, 'Color', [1 0 1], 'LineWidth', 1.5);
plot(tt, Es, 'Color', [1 0.9 0], 'LineWidth', 1.5);
plot(tt, Et, 'Color', [0 1 0.53], 'LineWidth', 1.5);
title('Ledger','Color','w'); xlabel('t'); dk(ax);
lgs(legend('kinetic','strain','total'));

try, set(fig,'InvertHardcopy','off'); print(fig,'fig16_15.png','-dpng','-r300'); catch, end
