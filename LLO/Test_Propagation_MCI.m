% Test_Propagation_MCI.m - Verify MCI perturbed two-body propagation
%
% This script tests the propagation of LLO orbits in the Moon-Centered Inertial (MCI) frame
% using the perturbed two-body dynamics (Moon gravity + Earth 3rd-body + J2).
%
% Generates a few LLO primaries, propagates them, and plots the results.

clear all
close all
clc

% Add paths
addpath(genpath('../../Population_generation'));
addpath(genpath('./'));

% Load defaults
cfg = LLO_defaults();
params = cfg.params;

% ODE options
ode_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% Time span for propagation (e.g., 1 day)
tspan = linspace(0, 24*3600, 1000); % 1 day, 1000 points

% Generate a few LLO primaries for testing
n_test = 3; % Number of orbits to test
families = {'Circular', 'Eccentric', 'Frozen'};
alt_range = [100, 200]; % km

fprintf('Generating test LLO orbits...\n');
test_states = cell(n_test, 1);
test_info = cell(n_test, 1);

for i = 1:n_test
    fam = families{mod(i-1, length(families)) + 1};
    fprintf('Generating %s orbit %d...\n', fam, i);
    
    % Generate one primary
    primaries_table = generate_LLO_primaries_MCI(1, fam, alt_range);
    
    % Extract state
    x0 = [primaries_table.x, primaries_table.y, primaries_table.z, ...
          primaries_table.vx, primaries_table.vy, primaries_table.vz]';
    
    test_states{i} = x0;
    test_info{i} = sprintf('%s (ID %d)', fam, primaries_table.ID);
end

% Propagate each orbit
fprintf('Propagating orbits...\n');
propagated_orbits = cell(n_test, 1);

for i = 1:n_test
    fprintf('Propagating %s...\n', test_info{i});
    [t, x] = ode113(@(t,x) two_body_perturbed_rhs(t, x, params), tspan, test_states{i}, ode_options);
    propagated_orbits{i} = x;
end

% Plot results
fprintf('Plotting results...\n');

% Plot 3D trajectories
figure('Position', [100, 100, 1200, 800]);
subplot(1, 2, 1);
hold on;
colors = {'b', 'r', 'g'};
for i = 1:n_test
    plot3(propagated_orbits{i}(:,1)/1000, propagated_orbits{i}(:,2)/1000, propagated_orbits{i}(:,3)/1000, ...
          colors{i}, 'LineWidth', 1.5);
end

% Plot Moon
[xm, ym, zm] = sphere(20);
surf(xm*cfg.params.R_moon_km, ym*cfg.params.R_moon_km, zm*cfg.params.R_moon_km, ...
     'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('LLO Trajectories in MCI Frame');
legend(test_info{:}, 'Location', 'best');
grid on;
axis equal;
view(3);

% Plot energy conservation
subplot(1, 2, 2);
hold on;
for i = 1:n_test
    r = vecnorm(propagated_orbits{i}(:,1:3), 2, 2);
    v = vecnorm(propagated_orbits{i}(:,4:6), 2, 2);
    E = v.^2 / 2 - params.GM_moon ./ r;
    plot(t/3600, E, colors{i}, 'LineWidth', 1.5);
end
xlabel('Time [hours]');
ylabel('Specific Energy [J/kg]');
title('Energy Conservation');
legend(test_info{:}, 'Location', 'best');
grid on;

% Save plot
saveas(gcf, 'Test_Propagation_MCI.png');
fprintf('Plot saved as Test_Propagation_MCI.png\n');

% Print final states
fprintf('\nFinal states after 1 day:\n');
for i = 1:n_test
    fprintf('%s:\n', test_info{i});
    fprintf('Position [km]: %.3f %.3f %.3f\n', propagated_orbits{i}(end,1:3)/1000);
    fprintf('Velocity [m/s]: %.3f %.3f %.3f\n', propagated_orbits{i}(end,4:6));
    r_final = norm(propagated_orbits{i}(end,1:3));
    fprintf('Altitude [km]: %.3f\n', (r_final - params.R_moon)/1000);
    fprintf('\n');
end

fprintf('Propagation test complete.\n');
