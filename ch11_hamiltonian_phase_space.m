% --- 1. The Hamiltonian of a pendulum: H = p^2 / (2 m l^2) - m g l cos(q) ---
m = 1; l = 1; g = 9.81;
H = @(q, p) p.^2 / (2 * m * l^2) - m * g * l * cos(q);

% Hamilton's equations: dq/dt = dH/dp,  dp/dt = -dH/dq
ham = @(t, s) [s(2) / (m * l^2); -m * g * l * sin(s(1))];

% --- 2. Phase-space landscape: every point (q, p) is a complete state of the system ---
q = linspace(-pi, pi, 600);                        % q is an angle, so phase space is a cylinder: show one turn
p = linspace(-9, 9, 400);
[Q, P] = meshgrid(q, p);

figure('Color', 'k'); hold on;
contour(Q, P, H(Q, P), 30, 'LineWidth', 0.6);
contour(Q, P, H(Q, P), [m * g * l, m * g * l], 'w', 'LineWidth', 1.5);   % the separatrix

% --- 3. Integrate a few trajectories: they stay on level sets of H ---
starts = [0.5 0; 1.5 0; 2.8 0; 0 7; 0 -7];         % rows are [q0, p0]
colors = {'c', 'm', 'y', 'g', [1 0.6 0]};
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);
worst = 0;
for i = 1:size(starts, 1)
    [~, S] = ode45(ham, linspace(0, 12, 3000), starts(i, :)', opts);
    Hs = H(S(:, 1), S(:, 2));
    worst = max(worst, max(abs(Hs - Hs(1))));
    qw = mod(S(:, 1) + pi, 2 * pi) - pi;               % wrap the angle into [-pi, pi)
    qw(find(abs(diff(qw)) > pi) + 1) = NaN;            % break the line at the wrap
    plot(qw, S(:, 2), 'Color', colors{i}, 'LineWidth', 2);
end
fprintf('largest energy drift along any trajectory: %.3e\n', worst);

set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
xlim([-pi pi]);
xlabel('position q (angle, wraps at +/- pi)', 'Color', 'w');
ylabel('momentum p', 'Color', 'w');
title('Hamiltonian phase space of a pendulum (white = separatrix)', 'Color', 'w');
