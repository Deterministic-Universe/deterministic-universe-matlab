% --- 1. The Lagrangian of a pendulum: L = K - U = (1/2) m l^2 th'^2 + m g l cos(th) ---
m = 1; l = 1; g = 9.81;
Lag = @(th, thd) 0.5 * m * l^2 * thd.^2 + m * g * l * cos(th);

% Euler-Lagrange: d/dt(dL/dth') - dL/dth = 0  ->  m l^2 th'' + m g l sin(th) = 0
pend = @(t, s) [s(2); -(g / l) * sin(s(1))];
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);

% --- 2. Solve the exact equation and the small-angle approximation ---
th0 = 2.5;                                         % a big swing (radians)
[t, S] = ode45(pend, linspace(0, 6, 6001), [th0; 0], opts);
small_angle = th0 * cos(sqrt(g / l) * t);

% --- 3. Check Euler-Lagrange numerically: d/dt(m l^2 th') + m g l sin(th) should be ~ 0 ---
p = m * l^2 * S(:, 2);                             % dL/dth' (the conjugate momentum)
residual = gradient(p, t) + m * g * l * sin(S(:, 1));
fprintf('max Euler-Lagrange residual: %.3e\n', max(abs(residual(3:end-2))));

% --- 4. Principle of least action: the true path makes S = integral of L dt stationary ---
T = 0.8;                                           % a short time window, where the stationary point is a minimum
[tt, tr] = ode45(pend, linspace(0, T, 4001), [0.6; 0], opts);
eps_values = linspace(-0.3, 0.3, 61);
action = zeros(size(eps_values));
for i = 1:numel(eps_values)                        % nudge the path, keeping both end points fixed
    e = eps_values(i);
    th  = tr(:, 1) + e * sin(pi * tt / T);
    thd = tr(:, 2) + e * (pi / T) * cos(pi * tt / T);
    action(i) = trapz(tt, Lag(th, thd));           % S = integral of L dt
end
[~, imin] = min(action);
fprintf('action is smallest at eps = %.2f\n', eps_values(imin));

% --- 5. Visual mapping ---
figure('Color', 'k');
subplot(1, 2, 1); hold on;
plot(t, S(:, 1), 'c', 'DisplayName', 'Euler-Lagrange (exact)');
plot(t, small_angle, 'm--', 'DisplayName', 'small-angle approximation');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Pendulum: exact vs small-angle', 'Color', 'w');
legend('Location', 'best');

subplot(1, 2, 2); hold on;
plot(eps_values, action, 'y');
yl = ylim; plot([0 0], yl, 'w');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
xlabel('\epsilon', 'Color', 'w');
title('Action S vs path perturbation: the true path is the minimum', 'Color', 'w');
