% Requires MATLAB R2016b or newer (local functions at the end of a script).

% --- 2. Initial Conditions: the figure-eight choreography (Chenciner & Montgomery, 2000) ---
G = 1;                                             % non-dimensionalized gravitational constant
masses = [1; 1; 1];                                % change these to leave the figure-eight (it becomes chaotic)
pos = [ 0.97000436, -0.24308753;                   % row i = position of body i
       -0.97000436,  0.24308753;
        0,           0         ];
v3 = [-0.93240737, -0.86473146];
vel = [-v3 / 2; -v3 / 2; v3] * sqrt(mean(masses)); % equal masses m: speeds scale as sqrt(m)

% Move to the center-of-mass frame so the system does not drift away
pos = pos - sum(masses .* pos, 1) / sum(masses);
vel = vel - sum(masses .* vel, 1) / sum(masses);

s0 = [reshape(pos', 6, 1); reshape(vel', 6, 1)];   % [x1;y1;x2;y2;x3;y3; vx1;vy1;...]
t_max = 10;
tspan = linspace(0, t_max, t_max * 500 + 1);
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, S] = ode45(@(t, s) nbody_rhs(t, s, masses, G), tspan, s0, opts);

% --- 3. Check the physics: total energy must stay constant ---
K = zeros(size(t));
for k = 1:3
    K = K + 0.5 * masses(k) * (S(:, 6 + 2*k - 1).^2 + S(:, 6 + 2*k).^2);
end
U = zeros(size(t));
pairs = [1 2; 1 3; 2 3];
for q = 1:3
    i = pairs(q, 1); j = pairs(q, 2);
    dist = sqrt((S(:, 2*i - 1) - S(:, 2*j - 1)).^2 + (S(:, 2*i) - S(:, 2*j)).^2);
    U = U - G * masses(i) * masses(j) ./ dist;
end
E = K + U;
fprintf('relative energy drift: %.3e\n', max(abs(E - E(1))) / abs(E(1)));

% --- 4. Generative Visual Mapping ---
colors = {'c', 'm', 'y'};
figure('Color', 'k');

% Visual Language 1: 'Shapes' (the orbits in ordinary space)
subplot(1, 2, 1); hold on;
for k = 1:3
    plot(S(:, 2*k - 1), S(:, 2*k), colors{k}, 'LineWidth', 0.8);
    plot(S(end, 2*k - 1), S(end, 2*k), 'o', 'MarkerFaceColor', colors{k}, 'MarkerEdgeColor', 'w', ...
         'MarkerSize', 6 + 3 * masses(k));
end
axis equal off;
title('Shapes (Orbits)', 'Color', 'w');

% Visual Language 2: 'Celestial' (gravitational potential at the final instant)
lim = 1.6;
[X, Y] = meshgrid(linspace(-lim, lim, 300));
Umap = zeros(size(X));
for k = 1:3
    Umap = Umap - G * masses(k) ./ sqrt((X - S(end, 2*k - 1)).^2 + (Y - S(end, 2*k)).^2 + 0.05^2);  % softened for display only
end
sorted = sort(Umap(:));
lo = sorted(round(0.02 * numel(sorted)));          % clip the deepest 2% so the colors stay readable
subplot(1, 2, 2);
contourf(X, Y, Umap, linspace(lo, max(Umap(:)), 50), 'LineStyle', 'none');
colormap(gca, hot);
axis equal off;
title('Celestial (Gravitational Potential Well)', 'Color', 'w');

% --- 1. The Physics Engine (local functions go at the end of a script) ---
function ds = nbody_rhs(~, s, masses, G)
    % s = [x1;y1;x2;y2;x3;y3; vx1;vy1;vx2;vy2;vx3;vy3]
    P = reshape(s(1:6), 2, 3)';                    % 3x2, row i = position of body i
    V = reshape(s(7:12), 2, 3)';                   % 3x2, row i = velocity of body i
    A = zeros(3, 2);
    for i = 1:3
        for j = 1:3
            if i ~= j
                d = P(j, :) - P(i, :);
                dist = max(norm(d), 1e-6);         % guard against a true singularity
                A(i, :) = A(i, :) + G * masses(j) * d / dist^3;   % a_i = sum_j G m_j (r_j - r_i)/|r_j - r_i|^3
            end
        end
    end
    ds = [reshape(V', 6, 1); reshape(A', 6, 1)];
end
