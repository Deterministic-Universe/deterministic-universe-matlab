% --- 1. A known path, so we can check the derivatives ---
t = linspace(0, 4, 2001);
r0 = [0; 0];
v0 = [12; 15];
a_true = [0; -9.81];
r = r0 + v0 .* t + 0.5 * a_true .* t.^2;           % 2 x N

% --- 2. Derivatives of space: v = dr/dt, a = dv/dt ---
v = [gradient(r(1, :), t); gradient(r(2, :), t)];  % first time derivative
a = [gradient(v(1, :), t); gradient(v(2, :), t)];  % second time derivative

% --- 3. Compare with the exact answers (skip the end points, where differences are one-sided) ---
v_exact = v0 + a_true .* t;
fprintf('max velocity error     : %.3e\n', max(max(abs(v(:, 2:end-1) - v_exact(:, 2:end-1)))));
fprintf('max acceleration error : %.3e\n', max(max(abs(a(:, 3:end-2) - a_true))));

% --- 4. Visual mapping ---
figure('Color', 'k');
subplot(1, 2, 1);
plot(t, v(1, :), 'c', t, v(2, :), 'm');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Velocity = dr/dt', 'Color', 'w'); legend('v_x', 'v_y');

subplot(1, 2, 2);
plot(t(3:end-2), a(1, 3:end-2), 'c', t(3:end-2), a(2, 3:end-2), 'm');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Acceleration = d^2r/dt^2', 'Color', 'w'); legend('a_x', 'a_y');
