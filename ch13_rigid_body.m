%% Chapter 13 -- Rigid Body Dynamics: Euler's Equations (Dzhanibekov effect)
% The Deterministic Universe -- Anthony Briggs
% Requires MATLAB R2016b or newer. Local functions are at the end of the file.
%
% Equations
%   I1 w1' = (I2 - I3) w2 w3
%   I2 w2' = (I3 - I1) w3 w1        (Euler, body frame)
%   I3 w3' = (I1 - I2) w1 w2
%   R' = R [w]x                     (orientation, body -> space)
%   E = 1/2 sum(Ik wk^2),  L_body = I w,  L_space = R * L_body
clear; clc; close all;

%% 1 Physics Engine
I = [1 2 3];                          % principal moments, I1 < I2 < I3
rhs = @(t, s) rigid_rhs(t, s, I);

%% 2 Initial Conditions
w0   = [0.05; 1.0; 0.05];            % spin about the INTERMEDIATE axis + tiny wobble
R0   = eye(3);
s0   = [w0; R0(:)];                  % state: [w1 w2 w3, 9 entries of R]
Tend = 60;
tspan = linspace(0, Tend, 2401);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, Y] = ode45(rhs, tspan, s0, opts);

w  = Y(:, 1:3)';                     % 3 x nt
nt = numel(t);
E  = 0.5 * sum(bsxfun(@times, I(:), w.^2), 1);
Lb = bsxfun(@times, I(:), w);        % body-frame angular momentum
Lmag = sqrt(sum(Lb.^2, 1));
Ls = zeros(3, nt);
for k = 1:nt
    R = reshape(Y(k, 4:12), 3, 3);
    Ls(:, k) = R * Lb(:, k);         % angular momentum in the fixed frame
end

%% 3 Generative Visual Mapping + ledger
fprintf('--- Checking your run (Chapter 13) ---\n');
fprintf('E(0)        = %.4f\n', E(1));
fprintf('|L|(0)      = %.4f\n', Lmag(1));
fprintf('max |E-E0|  = %.2e\n', max(abs(E - E(1))));
fprintf('max ||L|-L0| = %.2e\n', max(abs(Lmag - Lmag(1))));
fprintf('max drift of L_space = %.2e\n', max(max(abs(bsxfun(@minus, Ls, Ls(:,1))))));
flip = find(w(2,1:end-1) .* w(2,2:end) < 0);
fprintf('w2 changes sign at t = %s\n', mat2str(t(flip)', 4));
lambda = sqrt((I(2)-I(1)) * (I(3)-I(2)) / (I(1)*I(3))) * w0(2);
fprintf('linear growth rate about axis 2 (theory) = %.4f\n', lambda);

col = {[0 0.9 1], [1 0 1], [1 0.9 0]};
fig = figure('Color', 'k', 'Position', [100 100 1500 440]);

ax1 = subplot(1, 3, 1); hold on;
for k = 1:3, plot(t, w(k,:), 'Color', col{k}, 'LineWidth', 1.5); end
title('Body-frame angular velocity', 'Color', 'w'); xlabel('t');
lg = legend('\omega_1', '\omega_2', '\omega_3'); style_legend(lg); dark_axes(ax1);

ax2 = subplot(1, 3, 2); hold on;
for k = 1:3, plot(t, Ls(k,:), 'Color', col{k}, 'LineWidth', 1.5); end
title('Angular momentum in space (constant)', 'Color', 'w'); xlabel('t');
lg = legend('L_x', 'L_y', 'L_z'); style_legend(lg); dark_axes(ax2);

ax3 = subplot(1, 3, 3); hold on;
plot(t, E,    'Color', [0 1 0.53], 'LineWidth', 1.5);
plot(t, Lmag, 'Color', [1 0.53 0], 'LineWidth', 1.5);
ylim([0, 1.3 * max([E, Lmag])]);
title('Ledger: E and |L|', 'Color', 'w'); xlabel('t');
lg = legend('E', '|L|'); style_legend(lg); dark_axes(ax3);

set(fig, 'InvertHardcopy', 'off');
print(fig, 'fig16_13.png', '-dpng', '-r300');

%% ---------------- local functions ----------------
function ds = rigid_rhs(~, s, I)
    w = s(1:3);
    R = reshape(s(4:12), 3, 3);
    dw = [ (I(2) - I(3)) * w(2) * w(3) / I(1);
           (I(3) - I(1)) * w(3) * w(1) / I(2);
           (I(1) - I(2)) * w(1) * w(2) / I(3) ];
    S  = [  0    -w(3)  w(2);
            w(3)  0    -w(1);
           -w(2)  w(1)  0   ];
    dR = R * S;
    ds = [dw; dR(:)];
end

function dark_axes(ax)
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'Box', 'on');
end

function style_legend(lg)
    set(lg, 'Color', 'k', 'EdgeColor', 'w');
    try, set(lg, 'TextColor', 'w'); catch, end   % TextColor needs R2018a+
end
