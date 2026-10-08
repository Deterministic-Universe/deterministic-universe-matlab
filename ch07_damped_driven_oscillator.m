% --- 1. The Physics Engine: m x'' + c x' + k x = F0 cos(wd t) ---
m = 1; k = 4;                                      % natural frequency w0 = sqrt(k/m) = 2
osc = @(t, s, c, F0, wd) [s(2); (F0 * cos(wd * t) - c * s(2) - k * s(1)) / m];

% --- 2. Part A: free motion (F0 = 0) for four damping values ---
figure('Color', 'k');
subplot(1, 2, 1); hold on;
tspan = linspace(0, 12, 1200);
c_crit = 2 * sqrt(m * k);                          % critical damping
cs = [0, 0.4, c_crit, 8];
names = {'undamped', 'underdamped', 'critical', 'overdamped'};
colors = {'c', 'm', 'y', 'g'};
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
for i = 1:4
    [t, S] = ode45(@(t, s) osc(t, s, cs(i), 0, 0), tspan, [1; 0], opts);
    plot(t, S(:, 1), colors{i}, 'DisplayName', sprintf('%s (c = %.1f)', names{i}, cs(i)));
end
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Dissipation: the damping term c dx/dt', 'Color', 'w');
legend('Location', 'northeast');

% --- 3. Part B: resonance. Drive at many frequencies, measure the steady-state amplitude ---
c = 0.4; F0 = 1;
wd_curve = linspace(0.3, 4, 400);
exact = F0 ./ sqrt((k - m * wd_curve.^2).^2 + (c * wd_curve).^2);   % analytic steady-state amplitude
subplot(1, 2, 2); hold on;
plot(wd_curve, exact, 'w', 'DisplayName', 'exact formula');

wd_test = linspace(0.5, 3.5, 13);
amp = zeros(size(wd_test));
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
for i = 1:numel(wd_test)
    [t, S] = ode45(@(t, s) osc(t, s, c, F0, wd_test(i)), linspace(0, 120, 12001), [0; 0], opts);
    amp(i) = max(abs(S(t > 90, 1)));               % after the transient has died out
end
plot(wd_test, amp, 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'simulation');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Resonance: amplitude vs driving frequency', 'Color', 'w');
legend('Location', 'northeast');
exact_test = F0 ./ sqrt((k - m * wd_test.^2).^2 + (c * wd_test).^2);
fprintf('max |simulation - exact| = %.3e\n', max(abs(amp - exact_test)));
