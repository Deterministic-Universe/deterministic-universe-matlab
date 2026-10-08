% --- 1. The Physics Engine: F = -G M m / r^2 along r-hat  ->  a = -G M r / |r|^3 ---
% state s = [x; y; vx; vy] of a small body orbiting a heavy body fixed at the origin
GM = 1;                                            % dimensionless units (G = M = 1)
gravity = @(t, s) [s(3:4); -GM * s(1:2) / norm(s(1:2))^3];

% --- 2. Initial Conditions: start at r = 1 and vary only the launch speed ---
speeds = [1.0 0.8 1.2 1.5];
labels = {'Circular (v = 1.0)', 'Ellipse (v = 0.8)', 'Ellipse (v = 1.2)', 'Escape (v = 1.5 > sqrt(2))'};
colors = {'c', 'm', 'y', 'g'};
tspan = linspace(0, 14, 4000);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% --- 3. Generative Visual Mapping ---
figure('Color', 'k'); hold on;
for i = 1:numel(speeds)
    [~, S] = ode45(gravity, tspan, [1; 0; 0; speeds(i)], opts);
    plot(S(:, 1), S(:, 2), colors{i}, 'LineWidth', 1.2, 'DisplayName', labels{i});
end
plot(0, 0, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'HandleVisibility', 'off');
axis equal off; xlim([-3 3]); ylim([-3 3]);
title('Orbits under F = -GMm/r^2: launch speed sets the shape', 'Color', 'w');
legend('Location', 'southwest');
