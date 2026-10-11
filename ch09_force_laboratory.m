%% Chapter 9 -- The Force Laboratory (MATLAB)
% Self-contained: paste into the editor or the command window and run.
% No local functions, no other files needed. Needs MATLAB R2018a or newer.
% Every force law has the SAME signature:  F = law(t, r, v, p)
%   t = time, r = position [x; y], v = velocity [vx; vy], p = struct of constants
%   F = force on the particle, a 2-by-1 column vector.
% The engine adds up whichever laws are in the list:  m a = sum of the laws.
%   energy ledger: K = 1/2 m |v|^2,  U = -GM m / |r|,  flat K + U = conservative set
clear; clc; close all;

%% 1 Physics Engine -- the law library and the engine
p = struct('GM',1, 'm',1, 'qB',0.8, 'cd',0.5, 'c',0.15);

central  = @(t,r,v,p) -p.GM*p.m*r/max(norm(r),1e-6)^3;   % gravity-like pull toward the origin
magnetic = @(t,r,v,p)  p.qB*[v(2); -v(1)];                % magnetic force: at right angles to v
quadrag  = @(t,r,v,p) -p.cd*norm(v)*v;                    % quadratic drag
viscous  = @(t,r,v,p) -p.c*v;                             % viscous (linear) damping

% sum every active law:  each law returns a 2x1 column; sum them along dimension 2
total = @(laws,t,r,v,p) sum(cell2mat(cellfun(@(f) f(t,r,v,p), laws, 'UniformOutput', false)), 2);
% state s = [x; y; vx; vy];  the engine turns a list of laws into an ODE right-hand side
make  = @(laws) @(t,s) [s(3:4); total(laws,t,s(1:2),s(3:4),p)/p.m];
% ledger pieces for a state matrix with one row per time (columns x y vx vy)
Kof = @(S) 0.5*p.m*(S(:,3).^2 + S(:,4).^2);
Uof = @(S) -p.GM*p.m./max(hypot(S(:,1),S(:,2)),1e-6);

%% 2 Initial Conditions
s0    = [1; 0; 0; 0.9];                                   % r0 = (1, 0), v0 = (0, 0.9)
tspan = linspace(0, 40, 4001);
tight = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t, S] = ode45(make({central, magnetic}), tspan, s0, tight);

% extra runs for Q10: add a dissipative law. These orbits can fall toward the centre, so the
% central law is softened inside r = 0.2 (it never blows up) and each run is cut at r = 0.2.
softc = @(t,r,v,p) -p.GM*p.m*r/max(norm(r),0.2)^3;
loose = odeset('RelTol',1e-8,'AbsTol',1e-10);
[td, Sd] = ode45(make({softc, magnetic, quadrag}), tspan, s0, loose);
[tv, Sv] = ode45(make({softc, magnetic, viscous}), tspan, s0, loose);
kd = find(hypot(Sd(:,1),Sd(:,2)) < 0.2, 1);  if ~isempty(kd), td = td(1:kd-1); Sd = Sd(1:kd-1,:); end
kv = find(hypot(Sv(:,1),Sv(:,2)) < 0.2, 1);  if ~isempty(kv), tv = tv(1:kv-1); Sv = Sv(1:kv-1,:); end

K = Kof(S);  U = Uof(S);  E = K + U;
Ed = Kof(Sd) + Uof(Sd);  Ev = Kof(Sv) + Uof(Sv);
r  = hypot(S(:,1), S(:,2));

%% 3 Checks and figure
fprintf('--- Checking your run (Chapter 9) ---\n');
fprintf('E(0) = %.4f   max |E-E0| = %.2e\n', E(1), max(abs(E-E(1))));
fprintf('min K = %.4f   max U = %.4f   farthest r = %.3f   nearest r = %.3f\n', min(K), max(U), max(r), min(r));
Fc = central(0,s0(1:2),s0(3:4),p);  Fm = magnetic(0,s0(1:2),s0(3:4),p);
fprintf('At t = 0: central = (%.2f, %.2f)  magnetic = (%.2f, %.2f)  net a = (%.2f, %.2f)\n', ...
        Fc, Fm, (Fc+Fm)/p.m);
fprintf('dE/dt at t = 0:  quadratic drag = %.3f   viscous = %.3f\n', ...
        quadrag(0,s0(1:2),s0(3:4),p)'*s0(3:4), viscous(0,s0(1:2),s0(3:4),p)'*s0(3:4));
fprintf('E at the end:  conservative set %.3f   with quadratic drag %.3f   with viscous %.3f\n', E(end), Ed(end), Ev(end));

dk  = @(ax) set(ax,'Color','k','XColor','w','YColor','w','Box','on');
lgs = @(lg) set(lg,'Color','k','EdgeColor','w','TextColor','w');
fig = figure('Color','k','Position',[100 100 1300 500]);

ax = subplot(1,2,1); hold on;
plot(S(:,1), S(:,2), 'Color', [0 0.9 1], 'LineWidth', 1.4);
plot(S(1,1), S(1,2), 'o', 'Color', [0 1 0], 'MarkerFaceColor', [0 1 0]);       % start
plot(S(end,1), S(end,2), 'o', 'Color', [1 0 1], 'MarkerFaceColor', [1 0 1]);   % end
axis equal;
title('Trajectory under: central, magnetic','Color','w'); xlabel('x'); ylabel('y'); dk(ax);

ax = subplot(1,2,2); hold on;
plot(t, K, 'Color', [0 0.9 1],   'LineWidth', 1.5);
plot(t, U, 'Color', [1 0 1],     'LineWidth', 1.5);
plot(t, E, 'Color', [1 0.9 0],   'LineWidth', 1.5);
title('Energy ledger','Color','w'); xlabel('t'); ylabel('energy'); dk(ax);
lgs(legend('Kinetic','Potential (conservative laws)','K + U','Location','east'));

try, set(fig,'InvertHardcopy','off'); print(fig,'fig16_09.png','-dpng','-r300'); catch, end

% Optional Q10 picture: how the ledger bends when a dissipative law is added
fig2 = figure('Color','k','Position',[100 100 700 460]); ax = axes; hold on;
plot(t,  E,  'Color', [1 0.9 0],  'LineWidth', 1.5);
plot(td, Ed, 'Color', [1 0.3 0.3],'LineWidth', 1.5);
plot(tv, Ev, 'Color', [1 0.53 0], 'LineWidth', 1.5);
title('Dissipative laws bend the ledger down','Color','w'); xlabel('t'); ylabel('K + U'); dk(ax);
lgs(legend('conservative + magnetic','+ quadratic drag','+ viscous','Location','southwest'));
