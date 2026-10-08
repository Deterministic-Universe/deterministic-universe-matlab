% --- 1. The Physics Engine: r(t) = r0 + v0*t + (1/2)*a*t^2 ---
% Valid only for CONSTANT acceleration. r0, v0, a are 2x1 columns; t is a 1xN row.
position = @(t, r0, v0, a) r0 + v0 .* t + 0.5 * a .* t.^2;    % 2xN result

% --- 2. Initial Conditions ---
r0 = [0; 0];
speed = 20;
a = [0; -9.81];
t = linspace(0, 4.5, 450);

% --- 3. Generative Visual Mapping ---
figure('Color', 'k'); hold on;
for angle = 15:10:85
    v0 = speed * [cosd(angle); sind(angle)];
    r = position(t, r0, v0, a);
    r = r(:, r(2, :) >= 0);                        % stop at the ground
    plot(r(1, :), r(2, :), 'Color', [angle/90, 0.4, 1 - angle/90], 'LineWidth', 1.3);
end
axis equal off;
title('Parabolas from r(t) = r_0 + v_0 t + a t^2 / 2', 'Color', 'w');
