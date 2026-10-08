% Requires MATLAB R2016b or newer (local functions at the end of a script).

% --- 1. The Law Library: every law is a function handle returning a force on the particle ---
% Signature: law(t, r, v, p) -> 2x1 force.  r = position, v = velocity, p = parameters.
p = struct('m', 1, 'g', 9.81, 'k', 1, 'c', 0.15, 'cd', 0.1, 'F0', 0.5, 'wd', 1, 'GM', 1, 'qB', 0.8);
laws.gravity  = @(t, r, v, p) [0; -p.m * p.g];                        % F = m g
laws.spring   = @(t, r, v, p) -p.k * r;                               % Hooke: F = -k r
laws.damping  = @(t, r, v, p) -p.c * v;                               % F = -c v
laws.drag     = @(t, r, v, p) -p.cd * norm(v) * v;                    % F = -cd |v| v
laws.drive    = @(t, r, v, p) [p.F0 * cos(p.wd * t); 0];              % F0 cos(wd t)
laws.central  = @(t, r, v, p) -p.GM * p.m * r / max(norm(r), 1e-6)^3; % inverse-square attraction
laws.magnetic = @(t, r, v, p) p.qB * [v(2); -v(1)];                   % F = q v x B (does no work)

% --- 2. Choose the laws and the numbers: this is the part you play with ---
active = {'central', 'magnetic'};      % try {'spring','damping','drive'} or {'central','drag'}
r0 = [1; 0];
v0 = [0; 0.9];                         % with 'central', keep some sideways speed:
                                       % a particle that falls to r = 0 feels an infinite force
t_max = 40;

% --- 3. The Engine: sum whichever laws are active, then F = m a ---
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, S] = ode45(@(t, s) rhs(t, s, laws, active, p), linspace(0, t_max, 4000), [r0; v0], opts);

% --- 4. Energy ledger (potential energies of the conservative laws only) ---
pot.gravity = @(r, p) p.m * p.g * r(2);
pot.spring  = @(r, p) 0.5 * p.k * (r' * r);
pot.central = @(r, p) -p.GM * p.m / norm(r);
K = 0.5 * p.m * (S(:, 3).^2 + S(:, 4).^2);
U = zeros(size(K));
for i = 1:numel(active)
    if isfield(pot, active{i})         % drag, damping, drive and magnetic have no potential energy
        potential = pot.(active{i});
        for j = 1:numel(t)
            U(j) = U(j) + potential(S(j, 1:2)', p);
        end
    end
end
E = K + U;
fprintf('energy change over the run: %.3e\n', E(end) - E(1));

% --- 5. Visual mapping ---
figure('Color', 'k');
subplot(1, 2, 1); hold on;
plot(S(:, 1), S(:, 2), 'c', 'LineWidth', 1);
plot(S(1, 1), S(1, 2), 'go', 'MarkerFaceColor', 'g');
plot(S(end, 1), S(end, 2), 'mo', 'MarkerFaceColor', 'm');
axis equal off;
title(['Trajectory under: ' strjoin(active, ', ')], 'Color', 'w');

subplot(1, 2, 2); hold on;
plot(t, K, 'c', 'DisplayName', 'Kinetic');
plot(t, U, 'm', 'DisplayName', 'Potential (conservative laws)');
plot(t, E, 'y', 'LineWidth', 2, 'DisplayName', 'K + U');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title('Energy ledger (flat only if no damping, drag or drive)', 'Color', 'w');
legend('Location', 'best');

% --- The Engine's derivative function (local functions go at the end of a script) ---
function ds = rhs(t, s, laws, active, p)
    r = s(1:2);
    v = s(3:4);
    F = [0; 0];
    for i = 1:numel(active)
        law = laws.(active{i});
        F = F + law(t, r, v, p);
    end
    ds = [v; F / p.m];
end
