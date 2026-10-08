% --- 1. The Physics Engine: sum(F) = m * d2r/dt2  ->  a = F / m ---
% state s = [x; y; vx; vy];  d(s)/dt = [v; F/m]

% --- 2. Initial Conditions: the SAME force on four different masses ---
F = [0; 2];                                        % constant net force (N)
masses = [0.5 1.0 2.0 4.0];
s0 = [0; 0; 1; 0];                                 % start at the origin moving in +x
tspan = linspace(0, 4, 400);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
colors = {'c', 'm', 'y', 'g'};

% --- 3. Integrate and compare with the exact result r = r0 + v0 t + (F/2m) t^2 ---
figure('Color', 'k'); hold on;
for i = 1:numel(masses)
    m = masses(i);
    [t, S] = ode45(@(t, s) [s(3:4); F / m], tspan, s0, opts);
    exact_y = 0.5 * (F(2) / m) * t.^2;
    fprintf('m = %3.1f   max error vs exact = %.2e\n', m, max(abs(S(:, 2) - exact_y)));
    plot(S(:, 1), S(:, 2), colors{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('m = %g', m));
end
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Same force, different mass: heavier bodies curve less', 'Color', 'w');
legend('Location', 'best');
