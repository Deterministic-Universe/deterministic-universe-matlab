% --- 1. The Physics Engine: a small mass m orbiting a fixed heavy mass M ---
G = 1; M = 1; m = 1;
gravity = @(t, s) [s(3:4); -G * M * s(1:2) / norm(s(1:2))^3];

% --- 2. Integrate an elliptical orbit ---
tspan = linspace(0, 20, 4000);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, S] = ode45(gravity, tspan, [1; 0; 0; 0.8], opts);

% --- 3. Energy bookkeeping: E = K + U = (1/2) m v^2 - G M m / r ---
r = sqrt(S(:, 1).^2 + S(:, 2).^2);
v2 = S(:, 3).^2 + S(:, 4).^2;
K = 0.5 * m * v2;
U = -G * M * m ./ r;
E = K + U;
fprintf('E(0) = %.6f   relative drift = %.3e\n', E(1), max(abs(E - E(1))) / abs(E(1)));

% --- 4. Visual mapping: energy is traded, never created or destroyed ---
figure('Color', 'k'); hold on;
plot(t, K, 'c', 'DisplayName', 'Kinetic K');
plot(t, U, 'm', 'DisplayName', 'Potential U');
plot(t, E, 'y', 'LineWidth', 2, 'DisplayName', 'Total E = K + U');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
xlabel('time', 'Color', 'w');
title('Energy ledger of an orbit', 'Color', 'w');
legend('Location', 'best');
