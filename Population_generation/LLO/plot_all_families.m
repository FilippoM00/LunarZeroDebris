% PLOT_ALL_FAMILIES - Generate (and optionally propagate) LLOs from all families.
%
% Goal: Simulate how other pipelines use primaries generation — drive sampling
% entirely via LLO_defaults ranges (no hardcoded ranges here) and validate that
% the generated primaries respect those ranges. Propagation is optional.

clear; clc; close all;
% Ensure Population_generation root (and all subfolders) is on the path, independent of PoC_Analysis
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
popRoot = fileparts(thisDir); % ..\Population_generation
addpath(genpath(popRoot));

%% Define Physical Parameters (SI units, high precision)
cfg = LLO_defaults();
params = cfg.params;

% Simulation knobs (prefer cfg if present)
if isfield(cfg, 'n_primaries_per_family')
    n_orbits_per_family = cfg.n_primaries_per_family;
else
    n_orbits_per_family = 100; % default count per family
end

% Families to generate: prefer cfg.families_to_generate if provided
if isfield(cfg, 'families_to_generate')
    families_to_generate = cfg.families_to_generate;
else
    families_to_generate = {'Circular', 'Eccentric', 'Polar', 'Frozen', 'General', 'HighEccentric'};
end

% Propagation control (propagate by default, unless cfg.do_propagate = false)
if isfield(cfg, 'do_propagate')
    do_propagate = logical(cfg.do_propagate);
else
    do_propagate = true;
end
if ~do_propagate
    fprintf('Primaries generation validation mode (no propagation).\n');
end

%% Simulation Setup
% Propagation detail (only used if do_propagate = true)
n_periods = 3; % Longer propagation increases perturbation accumulation
points_per_period = 1000;

if do_propagate
    fprintf('⚠️  PERTURBATION ACCUMULATION WARNING:\n');
    fprintf('   - Propagating %d orbital periods\n', n_periods);
    fprintf('   - Perturbations (J2, Earth gravity) accumulate over time\n');
    fprintf('   - Frozen orbits may show secular variations in e and ω\n');
    fprintf('   - For TCA analysis, consider shorter propagation times\n\n');
end

% Print effective configuration summary (key ranges)
fprintf('Using LLO_defaults configuration for sampling:\n');
fprintf('  Circular: alt_range=[%g %g] km, inc=[%g %g] deg\n', cfg.alt_range(1), cfg.alt_range(2), cfg.circular_inc_range(1), cfg.circular_inc_range(2));
fprintf('  Eccentric: peri_alt=[%g %g] km, apo_alt=[%g %g] km, inc=[%g %g] deg\n', cfg.ecc_peri_alt_range(1), cfg.ecc_peri_alt_range(2), cfg.ecc_apo_alt_range(1), cfg.ecc_apo_alt_range(2), cfg.eccentric_inc_range(1), cfg.eccentric_inc_range(2));
fprintf('  Polar: alt_range=[%g %g] km, inc=[%g %g] deg\n', cfg.alt_range(1), cfg.alt_range(2), cfg.polar_inc_range(1), cfg.polar_inc_range(2));
fprintf('  Frozen families and weights:\n');
ff = fieldnames(cfg.frozen_families);
for ii=1:numel(ff)
    fam = cfg.frozen_families.(ff{ii});
    if isfield(fam, 'alt_range_km') && isfield(fam, 'ecc_range')
        fprintf('    %-12s: inc≈%.2f±%.2f deg, alt=[%g %g] km, e=[%g %g]\n', ff{ii}, fam.inc_deg, fam.inc_tolerance_deg, fam.alt_range_km(1), fam.alt_range_km(2), fam.ecc_range(1), fam.ecc_range(2));
    elseif isfield(fam, 'alt_range_km')
        fprintf('    %-12s: inc≈%.2f±%.2f deg, alt=[%g %g] km, e=(family default)\n', ff{ii}, fam.inc_deg, fam.inc_tolerance_deg, fam.alt_range_km(1), fam.alt_range_km(2));
    else
        fprintf('    %-12s: inc≈%.2f±%.2f deg, alt=(global), e=(family default)\n', ff{ii}, fam.inc_deg, fam.inc_tolerance_deg);
    end
end
fprintf('\n');

%% Generate Initial Conditions for All Families
fprintf('Generating primaries for all families...\n');
all_primaries = [];

for i = 1:length(families_to_generate)
    family_name = families_to_generate{i};
    fprintf('  - Generating %d %s LLOs...\n', n_orbits_per_family, family_name);
    % Pass [] for alt_range to let the generator use LLO_defaults internally
    primaries = generate_LLO_primaries_MCI(n_orbits_per_family, family_name, []);
    if isempty(all_primaries)
        all_primaries = primaries;
    else
        all_primaries = [all_primaries; primaries];
    end
end

n_total_orbits = height(all_primaries);
fprintf('Total orbits generated: %d\n', n_total_orbits);

%% Propagate and Store Orbits
% Validate sampled ranges by family
fprintf('\nValidating sampled ranges by family...\n');
fam_list = unique(all_primaries.Family);
for fi = 1:numel(fam_list)
    fname = fam_list{fi};
    sel = strcmp(all_primaries.Family, fname);
    alt = all_primaries.Altitude_km(sel);
    ecc = all_primaries.Eccentricity(sel);
    inc = all_primaries.Inclination_deg(sel);
    fprintf('  %-16s: N=%d, alt[min/max]=[%.1f %.1f] km, e[min/max]=[%.4f %.4f], i[min/max]=[%.2f %.2f] deg\n', ...
        fname, sum(sel), min(alt), max(alt), min(ecc), max(ecc), min(inc), max(inc));
end

% Prepare storage only if propagating
states_all = cell(n_total_orbits, 1);
orbital_elements_all = cell(n_total_orbits, 1); % Store orbital elements for stability analysis
if do_propagate
    fprintf('\nPropagating all orbits with perturbation monitoring...\n');
end

for k = 1:n_total_orbits
    % Get orbital elements for period calculation
    if strcmpi(all_primaries.Family{k}, 'Eccentric') || strcmpi(all_primaries.Family{k}, 'HighEccentric')
        periapsis_alt_km = all_primaries.Altitude_km(k);
        ecc = all_primaries.Eccentricity(k);
        r_p_m = params.R_moon + periapsis_alt_km * 1000;
        a_m = r_p_m / (1 - ecc);
    elseif strcmpi(all_primaries.Family{k}, 'General')
        % For general orbits, altitude is approximate - use as periapsis estimate
        alt_km = all_primaries.Altitude_km(k);
        ecc = all_primaries.Eccentricity(k);
        r_p_m = params.R_moon + alt_km * 1000;
        a_m = r_p_m / (1 - ecc);
    elseif startsWith(all_primaries.Family{k}, 'Frozen')
        % For frozen orbits, altitude is periapsis altitude
        periapsis_alt_km = all_primaries.Altitude_km(k);
        ecc = all_primaries.Eccentricity(k);
        r_p_m = params.R_moon + periapsis_alt_km * 1000;
        a_m = r_p_m / (1 - ecc);
    else % For Circular, Polar, altitude is based on semi-major axis
        alt_km = all_primaries.Altitude_km(k);
        a_m = params.R_moon + alt_km * 1000;
    end

    period_sec = 2 * pi * sqrt(a_m^3 / params.GM_moon);

    if do_propagate
        % Time span
        t_span = [0, n_periods * period_sec];
        t_vec = linspace(t_span(1), t_span(2), n_periods * points_per_period);

        % Initial state
        x0_MCI = [all_primaries.x(k); all_primaries.y(k); all_primaries.z(k); ...
                  all_primaries.vx(k); all_primaries.vy(k); all_primaries.vz(k)];

        % Propagate with higher precision options
        options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        [t_out, state_out] = ode45(@(t,s) two_body_perturbed_rhs(t, s, params), t_vec, x0_MCI, options);

        % Store trajectory in km
        states_all{k} = state_out(:, 1:3) / 1000;

        % Monitor orbital elements for stability analysis (especially for frozen orbits)
        if startsWith(all_primaries.Family{k}, 'Frozen') || strcmpi(all_primaries.Family{k}, 'General') || strcmpi(all_primaries.Family{k}, 'HighEccentric')
            orbital_elements = zeros(length(t_out), 6); % [a(km), e, i(deg), ω(deg), Ω(deg), ν(deg)]
            for j = 1:length(t_out)
                % Convert SI state to km-based units for canonical converter
                r_km   = state_out(j, 1:3)' / 1000;     % m -> km
                v_kms  = state_out(j, 4:6)' / 1000;     % m/s -> km/s
                mu_km3s2 = params.GM_moon / 1e9;        % m^3/s^2 -> km^3/s^2

                % Use canonical Utilities converter
                [a_km, e_val, i_rad, OM_rad, om_rad, th_rad] = car2kep(r_km, v_kms, mu_km3s2);

                % Store elements (angles in degrees)
                orbital_elements(j, :) = [a_km, e_val, i_rad*180/pi, om_rad*180/pi, OM_rad*180/pi, th_rad*180/pi];
            end
            orbital_elements_all{k} = orbital_elements;

            % Report perturbation effects for monitored orbits
            e_initial = orbital_elements(1, 2);
            e_final = orbital_elements(end, 2);
            e_drift = abs(e_final - e_initial);

            % Handle argument of perigee wrap-around (0° to 360°)
            omega_initial = orbital_elements(1, 4);
            omega_final = orbital_elements(end, 4);
            omega_diff = omega_final - omega_initial;

            % Handle wrap-around: if difference > 180°, subtract 360°; if < -180°, add 360°
            if omega_diff > 180
                omega_diff = omega_diff - 360;
            elseif omega_diff < -180
                omega_diff = omega_diff + 360;
            end
            omega_drift = abs(omega_diff);

            fprintf('  - Propagated orbit %d/%d (%s) - Perturbation effects:\n', k, n_total_orbits, all_primaries.Family{k});
            fprintf('    Eccentricity drift: %.2e (%.4f → %.4f)\n', e_drift, e_initial, e_final);
            fprintf('    Arg. perigee drift: %.2f° (%.2f° → %.2f°)\n', omega_drift, omega_initial, omega_final);
        else
            fprintf('  - Propagated orbit %d/%d (%s)\n', k, n_total_orbits, all_primaries.Family{k});
        end
    else
        % Not propagating: store just the initial point to enable a scatter plot
        states_all{k} = [all_primaries.x(k), all_primaries.y(k), all_primaries.z(k)]/1000;
    end
end

%% 3D Plot
figure('Name', 'Combined LLO Population - Perturbed 2-Body (MCI Frame)', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
hold on; axis equal; grid on;

% Plot Moon
[smx, smy, smz] = sphere(50);
surf(params.R_moon/1000*smx, params.R_moon/1000*smy, params.R_moon/1000*smz, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Define colors for each family
family_colors = struct('Circular', [0, 0.4470, 0.7410], ... % Blue
                       'Eccentric', [0.8500, 0.3250, 0.0980], ... % Orange
                       'Polar', [0.9290, 0.6940, 0.1250], ...    % Yellow
                       'Frozen', [0.4660, 0.6740, 0.1880], ...   % Green
                       'General', [0.4940, 0.1840, 0.5560], ...  % Purple
                       'HighEccentric', [0.6350, 0.0780, 0.1840], ... % Dark red
                       'Frozen_low_inc1', [0.3010, 0.7450, 0.9330], ... % Light blue
                       'Frozen_low_inc2', [0.6350, 0.0780, 0.1840], ... % Dark red
                       'Frozen_magic1', [0.4940, 0.1840, 0.5560], ...   % Purple
                       'Frozen_magic2', [0.9290, 0.6940, 0.1250], ...   % Yellow-orange
                       'Frozen_polar', [0.4660, 0.6740, 0.1880]);      % Green

% Plot orbits (as trajectories if propagated, else as points)
plot_handles = [];
for k = 1:n_total_orbits
    family = all_primaries.Family{k};
    color = family_colors.(family);
    if do_propagate
        h = plot3(states_all{k}(:,1), states_all{k}(:,2), states_all{k}(:,3), 'LineWidth', 1.5, 'Color', color);
    else
        h = plot3(states_all{k}(1), states_all{k}(2), states_all{k}(3), 'o', 'MarkerSize', 5, 'Color', color, 'MarkerFaceColor', color);
    end

    % Store one handle per family for the legend
    if ~isfield(plot_handles, family)
        plot_handles.(family) = h;
    end
end

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
if do_propagate
    title(sprintf('Combined LLO Population (%d orbits total)', n_total_orbits));
else
    title(sprintf('LLO Primaries (initial positions only) - %d samples', n_total_orbits));
end
view(45, 30);

% Create legend dynamically based on families present
unique_families = unique(all_primaries.Family);
legend_handles = [];
legend_names = {};
for i = 1:length(unique_families)
    family = unique_families{i};
    if isfield(plot_handles, family)
        legend_handles = [legend_handles, plot_handles.(family)];
        % Clean up family names for display
        if startsWith(family, 'Frozen_')
            display_name = strrep(family, 'Frozen_', 'Frozen ');
            display_name = strrep(display_name, '_', ' ');
        else
            display_name = family;
        end
        legend_names = [legend_names, display_name];
    end
end

legend(legend_handles, legend_names, 'Location', 'eastoutside');

fprintf('\nPlot generated.\n');
if do_propagate
    fprintf('⚠️  Remember: Perturbations accumulate over time.\n');
    fprintf('   For accurate TCA analysis, consider shorter propagation periods.\n');
end

%% Analyze Perturbation Effects Summary
if do_propagate
fprintf('\n📊 PERTURBATION EFFECTS ANALYSIS:\n');
fprintf('================================\n');

% Find frozen orbits and analyze their stability
frozen_indices = find(startsWith(all_primaries.Family, 'Frozen'));
general_indices = find(strcmpi(all_primaries.Family, 'General'));
highecc_indices = find(strcmpi(all_primaries.Family, 'HighEccentric'));

% Analyze frozen orbits
if ~isempty(frozen_indices)
    fprintf('Frozen Orbit Stability Analysis:\n');
    fprintf('-------------------------------\n');

    total_e_drift = 0;
    total_omega_drift = 0;
    n_frozen = 0;

    for idx = frozen_indices'
        if ~isempty(orbital_elements_all{idx})
            elements = orbital_elements_all{idx};
            e_drift = abs(elements(end, 2) - elements(1, 2));

            % Handle argument of perigee wrap-around
            omega_diff = elements(end, 4) - elements(1, 4);
            if omega_diff > 180
                omega_diff = omega_diff - 360;
            elseif omega_diff < -180
                omega_diff = omega_diff + 360;
            end
            omega_drift = abs(omega_diff);

            total_e_drift = total_e_drift + e_drift;
            total_omega_drift = total_omega_drift + omega_drift;
            n_frozen = n_frozen + 1;

            fprintf('  %s: e-drift=%.2e, ω-drift=%.2f°\n', ...
                   all_primaries.Family{idx}, e_drift, omega_drift);
        end
    end

    if n_frozen > 0
        avg_e_drift = total_e_drift / n_frozen;
        avg_omega_drift = total_omega_drift / n_frozen;

        fprintf('\nAverage perturbation effects over %d periods:\n', n_periods);
        fprintf('  Eccentricity drift: %.2e\n', avg_e_drift);
        fprintf('  Arg. perigee drift: %.2f°\n', avg_omega_drift);

        fprintf('\nStability Assessment:\n');
        if avg_e_drift < 1e-3 && avg_omega_drift < 10
            fprintf('  ✅ GOOD: Perturbations are minimal for TCA analysis\n');
        elseif avg_e_drift < 1e-2 && avg_omega_drift < 30
            fprintf('  ⚠️  MODERATE: Consider shorter propagation for precision\n');
        else
            fprintf('  ❌ SIGNIFICANT: Reduce propagation time for accurate TCA\n');
        end
    end
end

% Analyze general orbits
if ~isempty(general_indices)
    fprintf('\nGeneral Orbit Stability Analysis:\n');
    fprintf('--------------------------------\n');

    total_e_drift = 0;
    total_omega_drift = 0;
    n_general = 0;

    for idx = general_indices'
        if ~isempty(orbital_elements_all{idx})
            elements = orbital_elements_all{idx};
            e_drift = abs(elements(end, 2) - elements(1, 2));

            omega_diff = elements(end, 4) - elements(1, 4);
            if omega_diff > 180
                omega_diff = omega_diff - 360;
            elseif omega_diff < -180
                omega_diff = omega_diff + 360;
            end
            omega_drift = abs(omega_diff);

            total_e_drift = total_e_drift + e_drift;
            total_omega_drift = total_omega_drift + omega_drift;
            n_general = n_general + 1;

            fprintf('  General: e-drift=%.2e, ω-drift=%.2f°\n', e_drift, omega_drift);
        end
    end

    if n_general > 0
        fprintf('\nGeneral orbits average perturbation effects over %d periods:\n', n_periods);
        fprintf('  Eccentricity drift: %.2e\n', total_e_drift / n_general);
        fprintf('  Arg. perigee drift: %.2f°\n', total_omega_drift / n_general);
    end
end

% Analyze high eccentricity orbits
if ~isempty(highecc_indices)
    fprintf('\nHigh Eccentricity Orbit Stability Analysis:\n');
    fprintf('-----------------------------------------\n');

    total_e_drift = 0;
    total_omega_drift = 0;
    n_highecc = 0;

    for idx = highecc_indices'
        if ~isempty(orbital_elements_all{idx})
            elements = orbital_elements_all{idx};
            e_drift = abs(elements(end, 2) - elements(1, 2));

            omega_diff = elements(end, 4) - elements(1, 4);
            if omega_diff > 180
                omega_diff = omega_diff - 360;
            elseif omega_diff < -180
                omega_diff = omega_diff + 360;
            end
            omega_drift = abs(omega_diff);

            total_e_drift = total_e_drift + e_drift;
            total_omega_drift = total_omega_drift + omega_drift;
            n_highecc = n_highecc + 1;

            fprintf('  HighEcc: e-drift=%.2e, ω-drift=%.2f°\n', e_drift, omega_drift);
        end
    end

    if n_highecc > 0
        fprintf('\nHigh eccentricity orbits average perturbation effects over %d periods:\n', n_periods);
        fprintf('  Eccentricity drift: %.2e\n', total_e_drift / n_highecc);
        fprintf('  Arg. perigee drift: %.2f°\n', total_omega_drift / n_highecc);
    end
end

fprintf('\nRecommendations for TCA Analysis:\n');
fprintf('--------------------------------\n');
fprintf('1. For short-term TCA (< 1 orbit): Current setup is adequate\n');
fprintf('2. For medium-term TCA (1-10 orbits): Reduce n_periods to 1-2\n');
fprintf('3. For long-term TCA (> 10 orbits): Use simplified dynamics or\n');
fprintf('   implement orbit maintenance strategies\n');
fprintf('4. Monitor frozen orbit stability - perturbations accumulate!\n');
end
