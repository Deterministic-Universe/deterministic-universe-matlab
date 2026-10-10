%% Chapter 13 -- Rigid Body Dynamics (Euler's equations, Dzhanibekov flip)
% Self-contained: paste into the editor or the command window and run.
% No local functions, no other files needed. Needs MATLAB R2018a or newer.
%   I1 w1' = (I2-I3) w2 w3,  I2 w2' = (I3-I1) w3 w1,  I3 w3' = (I1-I2) w1 w2
%   R' = R [w]x,   E = 1/2 sum(Ik wk^2),   L_space = R (I w)
clear; clc; close all;

%% 1 Physics Engine
I   = [1 2 3];                                   % principal moments
col = @(A) A(:);                                 % turn a matrix into a column
sk  = @(w) [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];   % cross-product matrix [w]x
eul = @(s) [(I(2)-I(3))*s(2)*s(3)/I(1); ...
            (I(3)-I(1))*s(3)*s(1)/I(2); ...
            (I(1)-I(2))*s(1)*s(2)/I(3)];         % Euler's equations
rhs = @(t,s) [eul(s); col(reshape(s(4:12),3,3)*sk(s(1:3)))];

%% 2 Initial Conditions
w0 = [0.05; 1.0; 0.05];                          % spin about the intermediate axis + wobble
s0 = [w0; col(eye(3))];
tspan = linspace(0, 60, 2401);
[t, Y] = ode45(rhs, tspan, s0, odeset('RelTol',1e-10,'AbsTol',1e-12));

w  = Y(:,1:3)';                                  % 3 x nt
nt = numel(t);
E  = 0.5*sum(bsxfun(@times, I(:), w.^2), 1);
Lb = bsxfun(@times, I(:), w);
Lmag = sqrt(sum(Lb.^2, 1));
Ls = zeros(3, nt);
for k = 1:nt
    Ls(:,k) = reshape(Y(k,4:12),3,3) * Lb(:,k);  % angular momentum in the fixed frame
end

%% 3 Checks and figure
fprintf('--- Checking your run (Chapter 13) ---\n');
fprintf('E(0) = %.4f   |L|(0) = %.4f\n', E(1), Lmag(1));
fprintf('max |E-E0| = %.2e   max ||L|-L0| = %.2e\n', max(abs(E-E(1))), max(abs(Lmag-Lmag(1))));
fprintf('max drift of L in space = %.2e\n', max(max(abs(bsxfun(@minus, Ls, Ls(:,1))))));
flip = find(w(2,1:end-1).*w(2,2:end) < 0);
fprintf('w2 changes sign at t = %s\n', mat2str(t(flip)',4));
fprintf('linear growth rate about axis 2 = %.4f\n', sqrt((I(2)-I(1))*(I(3)-I(2))/(I(1)*I(3)))*w0(2));

c = {[0 0.9 1],[1 0 1],[1 0.9 0]};
dk  = @(ax) set(ax,'Color','k','XColor','w','YColor','w','Box','on');
lgs = @(lg) set(lg,'Color','k','EdgeColor','w','TextColor','w');
fig = figure('Color','k','Position',[100 100 1500 440]);

ax = subplot(1,3,1); hold on;
for k = 1:3, plot(t, w(k,:), 'Color', c{k}, 'LineWidth', 1.5); end
title('Body-frame angular velocity','Color','w'); xlabel('t'); dk(ax);
lgs(legend('\omega_1','\omega_2','\omega_3'));

ax = subplot(1,3,2); hold on;
for k = 1:3, plot(t, Ls(k,:), 'Color', c{k}, 'LineWidth', 1.5); end
title('Angular momentum in space (constant)','Color','w'); xlabel('t'); dk(ax);
lgs(legend('L_x','L_y','L_z'));

ax = subplot(1,3,3); hold on;
plot(t, E, 'Color', [0 1 0.53], 'LineWidth', 1.5);
plot(t, Lmag, 'Color', [1 0.53 0], 'LineWidth', 1.5);
ylim([0, 1.3*max([E, Lmag])]);
title('Ledger: E and |L|','Color','w'); xlabel('t'); dk(ax);
lgs(legend('E','|L|'));

try, set(fig,'InvertHardcopy','off'); print(fig,'fig16_13.png','-dpng','-r300'); catch, end
