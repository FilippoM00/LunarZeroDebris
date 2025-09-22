function [population, meta] = create_tca_populationLLO(total, opts)
%CREATE_TCA_POPULATIONLLO Entry point to generate a TCA population using LLO baseline.
%
%   [population, meta] = create_tca_populationLLO(total, opts)
%
%   Inputs
%     total (integer > 0): total number of close approaches to generate
%     opts  (struct, optional): configuration options
%
%   Outputs
%     population : struct with TCA population data
%     meta       : struct with metadata
%
%   This function adapts the TCA generation process to use the Verifiedv1.0
%   baseline with perturbed two-body dynamics in MCI frame and SI units.

% ======================================================================
% Input validation and initialization
% ======================================================================

% Start timing
tic;

% Validate total input
if nargin < 1 || isempty(total) || ~isreal(total) || ~isscalar(total) || total <= 0 || rem(total, 1) ~= 0
    error('Total must be a positive integer scalar');
end

% Handle optional opts argument
if nargin < 2 || isempty(opts)
    opts = struct();
end

% Validate opts structure
if ~isstruct(opts)
    error('opts must be a structure');
end

% Set default options
if ~isfield(opts, 'isotropic') || ~isstruct(opts.isotropic)
    opts.isotropic = struct();
end
if ~isfield(opts.isotropic, 'enable') || isempty(opts.isotropic.enable)
    % Match CR3BP behavior: isotropic azimuth by default
    opts.isotropic.enable = true;
end
if ~isfield(opts.isotropic, 'seed') || isempty(opts.isotropic.seed)
    opts.isotropic.seed = [];
end

% Add plotting and analysis options
if ~isfield(opts, 'plot') || ~isstruct(opts.plot)
    opts.plot = struct();
end
if ~isfield(opts.plot, 'enable') || isempty(opts.plot.enable)
    opts.plot.enable = false;
end
if ~isfield(opts.plot, 'saveDir') || isempty(opts.plot.saveDir)
    opts.plot.saveDir = fullfile(pwd, 'Media_LLO');
end
if ~isfield(opts.plot, 'format') || isempty(opts.plot.format)
    opts.plot.format = 'png';
end

% Ensure repository code and samplers are on MATLAB path (coherent with CR3BP flow)
here = fileparts(mfilename('fullpath'));
repo = fileparts(here); % repo root (parent of this file's directory)
addpath(genpath(repo));

% Load baseline defaults AFTER ensuring path, to avoid PoC shadowing
cfg = LLO_defaults();
params = cfg.params;

% Global RNG seeding (single point of control, coherent with CR3BP flow)
if isfield(opts,'seed') && ~isempty(opts.seed) && isreal(opts.seed) && ~isnan(opts.seed)
    rng(opts.seed);
end

% Validate configuration
if ~isfield(params, 'GM_moon') || isempty(params.GM_moon) || ~isreal(params.GM_moon) || params.GM_moon <= 0
    error('Invalid or missing GM_moon in configuration');
end
if ~isfield(params, 'R_moon') || isempty(params.R_moon) || ~isreal(params.R_moon) || params.R_moon <= 0
    error('Invalid or missing R_moon in configuration');
end
if ~isfield(cfg, 'alt_range') || isempty(cfg.alt_range) || ~isreal(cfg.alt_range) || any(cfg.alt_range <= 0)
    warning('Invalid altitude range in configuration. Using defaults.');
    cfg.alt_range = [50, 500]; % km
end

% Define families to generate (matching baseline)
families = {'Circular', 'Eccentric', 'Polar', 'Frozen', 'General', 'HighEccentric'};
n_families = numel(families);

% Generate primaries for each family (like plot_all_families.m)
primaries = struct();
total_primaries_generated = 0;

for fi = 1:n_families
    fam = families{fi};
    % Generate 100 primaries per family (same as original sampling)
    n_prim = 100;
    alt_range = cfg.alt_range; % Use default altitude range

    %fprintf('Generating %d %s LLO primaries...\n', n_prim, fam);
    if isfield(opts,'seed') && ~isempty(opts.seed) && isreal(opts.seed)
        primaries_table = generate_LLO_primaries_MCI(n_prim, fam, alt_range, struct('seed', opts.seed));
    else
        primaries_table = generate_LLO_primaries_MCI(n_prim, fam, alt_range);
    end

    % Validate generated primaries
    if isempty(primaries_table) || height(primaries_table) == 0
        error('No primaries generated for family %s', fam);
    end
    if height(primaries_table) ~= n_prim
        warning('Generated %d primaries for family %s, expected %d', height(primaries_table), fam, n_prim);
        n_prim = height(primaries_table);
    end

    % Convert to struct format similar to original
    primaries.(fam).states = [primaries_table.x, primaries_table.y, primaries_table.z, ...
                             primaries_table.vx, primaries_table.vy, primaries_table.vz]';
    primaries.(fam).count = n_prim;
    primaries.(fam).table = primaries_table; % Keep full table for metadata

    % Validate primary states
    if any(any(isnan(primaries.(fam).states))) || any(any(isinf(primaries.(fam).states)))
        error('Generated primaries for family %s contain NaN or Inf values', fam);
    end

    % Calculate proper orbital period for each orbit using vis-viva equation
    % Pre-allocate arrays for better performance
    r_mags = zeros(1, n_prim);
    v_mags = zeros(1, n_prim);
    energies = zeros(1, n_prim);
    periods = zeros(1, n_prim);

    for idx = 1:n_prim
        r_vec = [primaries_table.x(idx), primaries_table.y(idx), primaries_table.z(idx)];
        v_vec = [primaries_table.vx(idx), primaries_table.vy(idx), primaries_table.vz(idx)];
        r_mags(idx) = norm(r_vec);
        v_mags(idx) = norm(v_vec);
        energies(idx) = v_mags(idx)^2/2 - params.GM_moon/r_mags(idx);

        if energies(idx) >= 0
            % Unbound orbit, use circular approximation with current radius
            periods(idx) = 2*pi * sqrt(r_mags(idx)^3 / params.GM_moon);
        else
            % Bound orbit, calculate semi-major axis safely
            a = -params.GM_moon / (2 * energies(idx)); % semi-major axis
            if a > 0 && isreal(a) && ~isnan(a) && ~isinf(a)
                periods(idx) = 2*pi * sqrt(a^3 / params.GM_moon);
            else
                % Fallback for numerical issues
                periods(idx) = 2*pi * sqrt(r_mags(idx)^3 / params.GM_moon);
                warning('Invalid semi-major axis for orbit %d in family %s. Using circular approximation.', idx, fam);
            end
        end
    end
    primaries.(fam).T = periods;
    primaries.(fam).mu = repmat(params.GM_moon, 1, n_prim); % per-orbit mu

    total_primaries_generated = total_primaries_generated + n_prim;
end

fprintf('Total primaries generated: %d across %d families\n', total_primaries_generated, n_families);

% ======================================================================
% Step 2: Sample encounter parameters (same as original)
% ======================================================================

% Defaults for miss distance and time bounds
if ~isfield(opts,'miss') || ~isstruct(opts.miss), opts.miss = struct(); end
% Enforce a hard minimum floor of 1 mm (1e-6 km), even if user passes a smaller positive value
if ~isfield(opts.miss,'min_km') || isempty(opts.miss.min_km) || ~isreal(opts.miss.min_km) || isnan(opts.miss.min_km)
    opts.miss.min_km = 1e-6;
elseif opts.miss.min_km < 1e-6
    warning('Minimum miss distance below 1 mm provided (%.3e km). Clamping to 1e-6 km (1 mm).', opts.miss.min_km);
    opts.miss.min_km = 1e-6;
end
if ~isfield(opts.miss,'max_km') || isempty(opts.miss.max_km) || ~isreal(opts.miss.max_km) || isnan(opts.miss.max_km) || opts.miss.max_km <= opts.miss.min_km
    opts.miss.max_km = 3.0;
end
if ~isfield(opts,'time') || ~isstruct(opts.time), opts.time = struct(); end
if ~isfield(opts.time,'min_days') || isempty(opts.time.min_days) || ~isreal(opts.time.min_days) || isnan(opts.time.min_days) || opts.time.min_days < 0
    opts.time.min_days = 1.0;
end
if ~isfield(opts.time,'max_days') || isempty(opts.time.max_days) || ~isreal(opts.time.max_days) || isnan(opts.time.max_days) || opts.time.max_days <= opts.time.min_days
    opts.time.max_days = 4.0;
end

% Units conversion (SI units)
% TU_sec = 27.321661 * 86400; % seconds
% LU_m = 384400 * 1000; % meters

population.n_encounters = total;
meta.encounter = struct('by_family', struct());

% ======================================================================
% Step 3: Distribute encounters per family
% ======================================================================
n_per_family = floor(total / n_families);
remainder = rem(total, n_families);
counts_per_family = repmat(n_per_family, n_families, 1);
if remainder > 0
    counts_per_family(1:remainder) = counts_per_family(1:remainder) + 1;
end
meta.encounter_counts = struct('counts_per_family', counts_per_family, 'families', {families});

% ======================================================================
% Step 4: For each family, sample parameters and propagate
% ======================================================================
baseline_states = zeros(6, total);
baseline_mu = zeros(1, total);
baseline_family_idx = zeros(1, total);
baseline_primary_idx = zeros(1, total);
baseline_tau = zeros(1, total);

% Pre-allocate validation arrays
enc_types_valid = cell(total, 1);
vrel_valid = zeros(total, 1);
angles_valid = zeros(total, 1);
dmiss_valid = zeros(total, 1);
tca_valid = zeros(total, 1);

k = 0; % Initialize encounter counter

for fi = 1:n_families
    fam = families{fi};
    src = primaries.(fam);
    m = counts_per_family(fi);
    
    % Sample encounter parameters (same as original)
    encOpts = struct();
    if isfield(opts,'seed') && ~isempty(opts.seed) && isreal(opts.seed) && ~isnan(opts.seed)
        encOpts.seed = opts.seed;
    end
    try
        [types_f, vrel_f, angles_f, enc_meta_f] = sample_encounter_params_uniform(m, encOpts);
    catch ME
        warning(ME.identifier, 'Encounter parameter sampling failed for family %s: %s. Using defaults.', fam, ME.message);
        types_f = repmat({'head-on'}, m, 1);
        vrel_f = 0.1 * ones(m, 1); % 0.1 km/s
        angles_f = zeros(m, 1); % head-on
        enc_meta_f = struct();
    end
    meta.encounter.by_family.(fam) = enc_meta_f;
    
    % Sample TCA times
    timeOpts = struct('seed',123,'plot',false);
    try
        [tca_f, time_meta_f] = sample_tca_times_uniform(opts.time.min_days, opts.time.max_days, m, timeOpts);
        % Validate TCA times are within bounds and clamp if necessary
        tca_f = max(min(tca_f, opts.time.max_days), opts.time.min_days);
        % Ensure no NaN or Inf values
        if any(isnan(tca_f)) || any(isinf(tca_f))
            warning('TCA sampling produced invalid values for family %s. Using defaults.', fam);
            tca_f = opts.time.min_days + (opts.time.max_days - opts.time.min_days) * rand(m, 1);
        end
    catch ME
        warning(ME.identifier, 'TCA time sampling failed for family %s: %s. Using defaults.', fam, ME.message);
        tca_f = opts.time.min_days + (opts.time.max_days - opts.time.min_days) * rand(m, 1);
        time_meta_f = struct();
    end
    if ~isfield(meta,'time_by_family'), meta.time_by_family = struct(); end
    meta.time_by_family.(fam) = time_meta_f;
    
    % Sample miss distances
    try
        [dmiss_f, miss_meta_f] = sample_miss_distances_lognormal(opts.miss.min_km, opts.miss.max_km, m, ...
            struct('seed',42,'plot',false,'coverage',[0.02 0.98],'mode','resample'));
    catch ME
        warning(ME.identifier, 'Miss distance sampling failed for family %s: %s. Using defaults.', fam, ME.message);
        % Log-normal distribution approximation
        mu = log(opts.miss.min_km);
        sigma = (log(opts.miss.max_km) - log(opts.miss.min_km)) / 3; % 3-sigma range
        dmiss_f = exp(mu + sigma * randn(m, 1));
        dmiss_f = max(min(dmiss_f, opts.miss.max_km), opts.miss.min_km); % clamp to bounds
        miss_meta_f = struct();
    end
    if ~isfield(meta,'miss_by_family'), meta.miss_by_family = struct(); end
    meta.miss_by_family.(fam) = miss_meta_f;
    
    for j = 1:m
        idx = randi(src.count);
        % Validate index bounds
        if idx < 1 || idx > src.count || ~isreal(idx) || isnan(idx)
            warning('Invalid primary index %d for family %s. Using index 1.', idx, fam);
            idx = 1;
        end

        % Get initial state from generated primaries
        x0 = src.states(:, idx);
        mu_i = src.mu(idx); % per-orbit mu
        T = src.T(idx); % per-orbit period

        % PROPAGATE PRIMARY FOR RANDOM TIME (TCA-like conditions)
        % This simulates the primary being at a random point in its orbit during encounter
        tau = rand() * T; % Random time between 0 and T (like TCA timing)

        % Validate propagation time
        if tau <= 0 || tau > T || ~isreal(tau) || isnan(tau) || isinf(tau)
            warning('Invalid propagation time for family %s, orbit %d. Using T/2.', fam, idx);
            tau = T/2;
        end

        % Validate initial state
        if any(isnan(x0)) || any(isinf(x0))
            warning('Initial state contains NaN or Inf for family %s, orbit %d. Skipping.', fam, idx);
            continue; % Skip this orbit, don't increment k
        end

        k = k + 1; % Only increment k when we have valid data to store

        % Update ODE options with current tau
        ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12, 'Stats','off', 'OutputFcn',[], ...
                           'MaxStep', tau/10, 'InitialStep', tau/100);

        try
            [~, X] = ode113(@(t,x) two_body_perturbed_rhs(t,x,params), [0, tau], x0(:)', ode_options);
            x_tau = X(end, :)';
        catch ME
            warning(ME.identifier, 'Primary propagation failed for family %s, orbit %d: %s. Attempting fallback.', fam, idx, ME.message);
            % Fallback: try with smaller time step
            try
                dt_small = min(tau/10, 60); % 60 seconds max step
                t_span = 0:dt_small:tau;
                if abs(t_span(end) - tau) > 1e-6
                    t_span = [t_span, tau];
                end
                [~, X] = ode113(@(t,x) two_body_perturbed_rhs(t,x,params), t_span, x0(:)', ode_options);
                x_tau = X(end, :)';
            catch ME2
                warning(ME2.identifier, 'Fallback propagation also failed: %s. Skipping this encounter.', ME2.message);
                % Decrement k since we couldn't generate this encounter
                k = k - 1;
                continue; % Skip this encounter entirely
            end
        end

        % Validate propagated primary state
        if any(isnan(x_tau)) || any(isinf(x_tau))
            warning('Propagated primary state contains NaN or Inf for family %s, orbit %d. Using initial state.', fam, idx);
            x_tau = x0;
        elseif norm(x_tau(1:3)) < params.R_moon
            warning('Propagated primary intersects Moon surface for family %s, orbit %d. Using initial state.', fam, idx);
            x_tau = x0;
        end

        % Store the propagated primary state as baseline
        baseline_states(:, k) = x_tau;
        baseline_mu(k) = mu_i;
        baseline_family_idx(k) = fi;
        baseline_primary_idx(k) = idx;
        baseline_tau(k) = tau; % Store the random propagation time

        % Validate baseline state
        r_baseline = norm(x_tau(1:3));
        if r_baseline < params.R_moon
            warning('Baseline state intersects Moon surface for encounter %d (r=%.1f km). This should not happen.', k, r_baseline/1000);
        end

        % Attach stratified parameters
        enc_types_valid{k} = types_f{j};
        vrel_valid(k) = vrel_f(j);
        angles_valid(k) = angles_f(j);
        dmiss_valid(k) = dmiss_f(j);
        tca_valid(k) = tca_f(j);
    end
end

% Set actual number of encounters generated
N = k;

% Check if we generated any encounters
if N == 0
    error('No valid encounters could be generated. Check input parameters and primary generation.');
end

% Attach baseline outputs
population.baseline_states = baseline_states(:, 1:N);
population.baseline_mu = baseline_mu(1:N);
population.baseline_family_idx = baseline_family_idx(1:N);
population.baseline_primary_idx = baseline_primary_idx(1:N);
population.baseline_tau = baseline_tau(1:N);
population.encounter_types = enc_types_valid(1:N);
population.vrel_kms = vrel_valid(1:N);
population.angles_deg = angles_valid(1:N);
population.dmiss_km = dmiss_valid(1:N);
population.tca_days = tca_valid(1:N);

% Add missing population fields to match original
population.primaries = primaries;
population.total_requested = total;
population.total_generated = N; % Actual number generated (may be less due to skips)
population.families = families;
secondary_states = zeros(6, N);
r_rel_all = zeros(3, N);
v_rel_all = zeros(3, N);

% Prepare RNG for isotropic azimuths once (coherent with CR3BP)
if isfield(opts, 'isotropic') && isfield(opts.isotropic, 'seed') && ~isempty(opts.isotropic.seed)
    rng(opts.isotropic.seed);
end

for k2 = 1:N
    x_p = baseline_states(:, k2);
    r_p = x_p(1:3);
    v_p = x_p(4:6);
    
    % Local encounter frame in MCI
    v_norm = norm(v_p);
    if v_norm < 1e-12
        % Fallback frame for near-zero velocity (e.g., circular orbits at apoapsis)
        e_r = r_p / max(norm(r_p), eps);
        % Use a more robust fallback that doesn't assume coordinate system orientation
        e_t = cross([0;0;1], e_r); % Cross with z-axis
        if norm(e_t) < 1e-12
            e_t = cross([1;0;0], e_r); % If z-axis aligned, use x-axis
        end
        e_t = e_t / max(norm(e_t), eps);
        e_n = cross(e_t, e_r); e_n = e_n / max(norm(e_n), eps);
        e_l = cross(e_n, e_t); e_l = e_l / max(norm(e_l), eps);
    else
        e_t = v_p / v_norm;
        e_r = r_p / max(norm(r_p), eps);
        e_n = cross(e_t, e_r); e_n = e_n / max(norm(e_n), eps);
        e_l = cross(e_n, e_t); e_l = e_l / max(norm(e_l), eps);
    end
    
    % v_rel direction
    theta = angles_valid(k2);
    u_v = cosd(theta)*e_t + sind(theta)*e_l;
    u_v = u_v / max(norm(u_v), eps);
    
    % Magnitudes in SI units
    vrel_ms = vrel_valid(k2) * 1000; % km/s to m/s
    dmiss_m = dmiss_valid(k2) * 1000; % km to m
    
    v_rel_vec = vrel_ms * u_v;
    
    % Perpendicular plane basis
    e_b1 = cross(u_v, e_n);
    if norm(e_b1) < 1e-12
        e_b1 = cross(u_v, e_t);
        if norm(e_b1) < 1e-12
            e_b1 = cross(u_v, e_r); % Final fallback
        end
    end
    if norm(e_b1) < 1e-12
        warning('Unable to construct perpendicular basis for encounter %d. Using default.', k2);
        e_b1 = [1;0;0]; % Default fallback
    end
    e_b1 = e_b1 / max(norm(e_b1), eps);
    e_b2 = cross(u_v, e_b1);
    e_b2 = e_b2 / max(norm(e_b2), eps);
    
    if opts.isotropic.enable
        phi = 2*pi*rand();
        b_hat = cos(phi)*e_b1 + sin(phi)*e_b2;
    else
        b_hat = e_b1;
    end
    r_rel_vec = dmiss_m * b_hat;
    
    % Secondary state
    x_s = [r_p + r_rel_vec; v_p + v_rel_vec];
    secondary_states(:, k2) = x_s;
    r_rel_all(:, k2) = r_rel_vec;
    v_rel_all(:, k2) = v_rel_vec;
    
    % Validate secondary state
    r_secondary = norm(x_s(1:3));
    if r_secondary < params.R_moon
        % Preserve TCA constraints: rotate r_rel within the plane ⟂ u_v towards outward radial
        warning('Secondary inside Moon for encounter %d (r=%.1f km). Re-orienting r_rel outward.', k2, r_secondary/1000);
        r_hat = r_p / max(norm(r_p), eps);
        b_out = r_hat - dot(r_hat, u_v)*u_v; % project radial onto plane ⟂ u_v
        if norm(b_out) > 1e-12
            b_hat = b_out / norm(b_out);
            r_rel_vec = dmiss_m * b_hat;
            x_s = [r_p + r_rel_vec; v_p + v_rel_vec];
            secondary_states(:, k2) = x_s;
            r_rel_all(:, k2) = r_rel_vec;
        end
        % If still inside, keep state but issue warning; user may resample
        if norm(x_s(1:3)) < params.R_moon
            warning('Re-oriented secondary still inside Moon for encounter %d. Consider resampling azimuth.', k2);
        end
    end
end

% Attach secondary states to population
population.secondary_states = secondary_states;
population.r_rel = r_rel_all;
population.v_rel = v_rel_all;

% ======================================================================
% Enhanced Characterization and Plotting (Optional)
% ======================================================================

if opts.plot.enable
    fprintf('Generating enhanced characterization plots...\n');
    
    % Create media directory if it doesn't exist
    if ~exist(opts.plot.saveDir, 'dir')
        mkdir(opts.plot.saveDir);
        fprintf('Created media directory: %s\n', opts.plot.saveDir);
    end
    
    % Generate comprehensive plots
    try
        % 1. Miss Distance Distribution Plot
        plot_miss_distance_analysis(population.dmiss_km, opts);
        
        % 2. TCA Times Distribution Plot  
        plot_tca_times_analysis(population.tca_days, opts);
        
        % 3. Relative Velocity Distribution Plot
        plot_relative_velocity_analysis(population.vrel_kms, opts);
        
        % 4. Encounter Angles Distribution Plot
        plot_encounter_angles_analysis(population.angles_deg, opts);
        
        % 5. Encounter Types Distribution Plot
        plot_encounter_types_analysis(population.encounter_types, opts);
        
        % 6. Family-wise Statistics Plot
        plot_family_statistics(population, families, opts);
        
        % 7. Orbital Elements Analysis
        plot_orbital_elements_analysis(population.baseline_states, population.baseline_mu, opts);
        
        fprintf('All characterization plots generated successfully.\n');
        
    catch ME
        warning(ME.identifier, 'Plotting failed: %s. Continuing with analysis.', ME.message);
    end
    
    % Generate statistical summary
    try
        stats_summary = generate_population_statistics(population, families);
        meta.statistics = stats_summary;
        fprintf('Statistical analysis completed.\n');
        
        % Print summary to console
        print_population_summary(stats_summary);
        
    catch ME
        warning(ME.identifier, 'Statistical analysis failed: %s', ME.message);
    end
end

% ======================================================================
% Final validation and summary
% ======================================================================

% Validate final population
if size(population.baseline_states, 2) ~= N
    error('Population size mismatch: expected %d, got %d', N, size(population.baseline_states, 2));
end

if size(population.secondary_states, 2) ~= N
    error('Secondary states size mismatch: expected %d, got %d', N, size(population.secondary_states, 2));
end

% Summary statistics
meta.summary = struct();
meta.summary.total_requested = total;
meta.summary.total_generated = N;
meta.summary.encounters_skipped = total - N;
meta.summary.families_processed = n_families;
meta.summary.primary_orbits_per_family = n_prim;
meta.summary.propagation_success_rate = (N / total) * 100;
meta.summary.average_miss_distance_km = mean(dmiss_valid(1:N));
meta.summary.average_vrel_kms = mean(vrel_valid(1:N));
meta.summary.generation_timestamp = datestr(now);

% ======================================================================
% Quality validation (enhanced from original create_tca_population.m)
% ======================================================================

fprintf('Performing quality validation...\n');

% Orthogonality and magnitude checks (like original)
orthogonality_errors = zeros(1, N);
miss_distance_errors = zeros(1, N);
vrel_errors = zeros(1, N);

for q = 1:N
    % Check orthogonality: r_rel · v_rel should be ≈ 0
    orthogonality_errors(q) = abs(dot(r_rel_all(:, q), v_rel_all(:, q)));
    
    % Check miss distance accuracy
    r_rel_magnitude = norm(r_rel_all(:, q));
    miss_distance_errors(q) = abs(r_rel_magnitude/1000 - dmiss_valid(q));
    
    % Check relative velocity accuracy
    v_rel_magnitude = norm(v_rel_all(:, q));
    vrel_errors(q) = abs(v_rel_magnitude/1000 - vrel_valid(q));
end

meta.validation = struct();
meta.validation.orthogonality_max = max(orthogonality_errors);
meta.validation.orthogonality_mean = mean(orthogonality_errors);
meta.validation.orthogonality_std = std(orthogonality_errors);
meta.validation.miss_distance_km_max_error = max(miss_distance_errors);
meta.validation.miss_distance_km_mean_error = mean(miss_distance_errors);
meta.validation.miss_distance_km_std_error = std(miss_distance_errors);
meta.validation.vrel_kms_max_error = max(vrel_errors);
meta.validation.vrel_kms_mean_error = mean(vrel_errors);
meta.validation.vrel_kms_std_error = std(vrel_errors);
meta.validation.total_encounters_validated = N;
meta.validation.validation_timestamp = datestr(now);

fprintf('Validation complete:\n');
fprintf('  Max orthogonality error: %.2e\n', meta.validation.orthogonality_max);
fprintf('  Max miss distance error: %.2e km\n', meta.validation.miss_distance_km_max_error);
fprintf('  Max relative velocity error: %.2e km/s\n', meta.validation.vrel_kms_max_error);

% ======================================================================
% Enhanced metadata structure (matching original create_tca_population.m)
% ======================================================================

meta.stage = 'complete';
meta.version = '2.0';
meta.framework = 'Verifiedv1.0';
computation_time = toc;
meta.computation_time = computation_time;
meta.summary.computation_time_seconds = computation_time;

% Check actual family counts
if isstruct(primaries)
    meta.counts_actual = struct();
    family_names = fieldnames(primaries);
    for fi = 1:numel(family_names)
        fam = family_names{fi};
        if isfield(primaries.(fam), 'count')
            meta.counts_actual.(fam) = primaries.(fam).count;
        else
            meta.counts_actual.(fam) = 'count field missing';
        end
    end
end
meta.timestamp = datestr(now);
meta.matlab_version = version();
meta.notes = sprintf('Generated %d TCA cases across %d LLO families with on-the-fly primary generation and random TCA timing. Enhanced validation and characterization.', N, n_families);

end

% ======================================================================
% Helper Functions for Enhanced Characterization and Plotting
% ======================================================================

function plot_miss_distance_analysis(dmiss_km, opts)
% Plot miss distance distribution analysis
figure('Position', [100, 100, 800, 600]);

% Log-normal histogram
subplot(2, 2, 1);
histogram(log10(dmiss_km), 40, 'Normalization', 'pdf');
xlabel('log_{10}(Miss Distance) [km]');
ylabel('Probability Density');
title('Miss Distance Distribution (Log Scale)');
grid on;

% Linear histogram  
subplot(2, 2, 2);
histogram(dmiss_km, 40, 'Normalization', 'pdf');
xlabel('Miss Distance [km]');
ylabel('Probability Density');
title('Miss Distance Distribution (Linear Scale)');
grid on;

% Cumulative distribution
subplot(2, 2, 3);
[f, x] = ecdf(dmiss_km);
plot(x, f, 'LineWidth', 2);
xlabel('Miss Distance [km]');
ylabel('Cumulative Probability');
title('Miss Distance CDF');
grid on;

% Statistics text
subplot(2, 2, 4);
axis off;
stats_text = sprintf(['Miss Distance Statistics:\n' ...
                     'Min: %.2e km\n' ...
                     'Max: %.2e km\n' ...
                     'Mean: %.2e km\n' ...
                     'Median: %.2e km\n' ...
                     'Std: %.2e km\n' ...
                     'Total Samples: %d'], ...
                     min(dmiss_km), max(dmiss_km), mean(dmiss_km), ...
                     median(dmiss_km), std(dmiss_km), length(dmiss_km));
text(0.1, 0.9, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'miss_distance_analysis');
end

function plot_tca_times_analysis(tca_days, opts)
% Plot TCA times distribution analysis
figure('Position', [100, 100, 800, 600]);

% Histogram
subplot(2, 2, 1);
histogram(tca_days, 30, 'Normalization', 'pdf');
xlabel('TCA Time [days]');
ylabel('Probability Density');
title('TCA Times Distribution');
grid on;

% Cumulative distribution
subplot(2, 2, 2);
[f, x] = ecdf(tca_days);
plot(x, f, 'LineWidth', 2);
xlabel('TCA Time [days]');
ylabel('Cumulative Probability');
title('TCA Times CDF');
grid on;

% Time series plot (sorted)
subplot(2, 2, 3);
plot(sort(tca_days), 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('TCA Time [days]');
title('TCA Times (Sorted)');
grid on;

% Statistics text
subplot(2, 2, 4);
axis off;
stats_text = sprintf(['TCA Times Statistics:\n' ...
                     'Min: %.3f days\n' ...
                     'Max: %.3f days\n' ...
                     'Mean: %.3f days\n' ...
                     'Median: %.3f days\n' ...
                     'Std: %.3f days\n' ...
                     'Total Samples: %d'], ...
                     min(tca_days), max(tca_days), mean(tca_days), ...
                     median(tca_days), std(tca_days), length(tca_days));
text(0.1, 0.9, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'tca_times_analysis');
end

function plot_relative_velocity_analysis(vrel_kms, opts)
% Plot relative velocity distribution analysis
figure('Position', [100, 100, 800, 600]);

% Linear histogram
subplot(2, 2, 1);
histogram(vrel_kms, 40, 'Normalization', 'pdf');
xlabel('Relative Velocity [km/s]');
ylabel('Probability Density');
title('Relative Velocity Distribution');
grid on;

% Log density plot
subplot(2, 2, 2);
log_vrel = log10(max(vrel_kms, 1e-10));
histogram(log_vrel, 40, 'Normalization', 'pdf');
xlabel('log_{10}(Relative Velocity) [km/s]');
ylabel('Probability Density');
title('Relative Velocity Distribution (Log Scale)');
grid on;

% Cumulative distribution
subplot(2, 2, 3);
[f, x] = ecdf(vrel_kms);
plot(x, f, 'LineWidth', 2);
xlabel('Relative Velocity [km/s]');
ylabel('Cumulative Probability');
title('Relative Velocity CDF');
grid on;

% Statistics text
subplot(2, 2, 4);
axis off;
stats_text = sprintf(['Relative Velocity Statistics:\n' ...
                     'Min: %.4f km/s\n' ...
                     'Max: %.4f km/s\n' ...
                     'Mean: %.4f km/s\n' ...
                     'Median: %.4f km/s\n' ...
                     'Std: %.4f km/s\n' ...
                     'Total Samples: %d'], ...
                     min(vrel_kms), max(vrel_kms), mean(vrel_kms), ...
                     median(vrel_kms), std(vrel_kms), length(vrel_kms));
text(0.1, 0.9, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'relative_velocity_analysis');
end

function plot_encounter_angles_analysis(angles_deg, opts)
% Plot encounter angles distribution analysis
figure('Position', [100, 100, 800, 600]);

% Angular histogram
subplot(2, 2, 1);
histogram(angles_deg, 36, 'Normalization', 'pdf');
xlabel('Encounter Angle [deg]');
ylabel('Probability Density');
title('Encounter Angles Distribution');
grid on;

% Polar plot
subplot(2, 2, 2);
polarhistogram(deg2rad(angles_deg), 36, 'Normalization', 'pdf');
title('Encounter Angles (Polar)');

% Cumulative distribution
subplot(2, 2, 3);
[f, x] = ecdf(angles_deg);
plot(x, f, 'LineWidth', 2);
xlabel('Encounter Angle [deg]');
ylabel('Cumulative Probability');
title('Encounter Angles CDF');
grid on;

% Statistics text
subplot(2, 2, 4);
axis off;
stats_text = sprintf(['Encounter Angles Statistics:\n' ...
                     'Min: %.1f°\n' ...
                     'Max: %.1f°\n' ...
                     'Mean: %.1f°\n' ...
                     'Median: %.1f°\n' ...
                     'Std: %.1f°\n' ...
                     'Total Samples: %d'], ...
                     min(angles_deg), max(angles_deg), mean(angles_deg), ...
                     median(angles_deg), std(angles_deg), length(angles_deg));
text(0.1, 0.9, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'encounter_angles_analysis');
end

function plot_encounter_types_analysis(encounter_types, opts)
% Plot encounter types distribution analysis
figure('Position', [100, 100, 1200, 800]);

% Count unique types
[unique_types, ~, ic] = unique(encounter_types);
counts = accumarray(ic, 1);

% Bar chart (uniform order)
subplot(2, 3, 1);
bar(counts);
set(gca, 'XTickLabel', unique_types);
xlabel('Encounter Type');
ylabel('Count');
title('Encounter Types (Uniform Order)');
grid on;

% Pie chart
subplot(2, 3, 2);
pie(counts, unique_types);
title('Encounter Types (Pie Chart)');

% Pareto chart (sorted by frequency - most common first)
subplot(2, 3, 3);
[sorted_counts, sort_idx] = sort(counts, 'descend');
sorted_types = unique_types(sort_idx);
bar(sorted_counts);
set(gca, 'XTickLabel', sorted_types);
xlabel('Encounter Type');
ylabel('Count');
title('Encounter Types (Sorted by Frequency)');
grid on;

% Uniform distribution chart (original order)
subplot(2, 3, 4);
bar(counts);
set(gca, 'XTickLabel', unique_types);
xlabel('Encounter Type');
ylabel('Count');
title('Encounter Types (Uniform Order)');
grid on;

% Statistics summary
subplot(2, 3, 5);
axis off;
stats_text = sprintf(['Encounter Types Summary:\n' ...
                     'Total Types: %d\n' ...
                     'Total Samples: %d\n\n' ...
                     'Top Types:\n'], ...
                     length(unique_types), length(encounter_types));
for i = 1:min(5, length(sorted_types))
    stats_text = sprintf('%s%s: %d (%.1f%%)\n', ...
                        stats_text, sorted_types{i}, sorted_counts(i), ...
                        100*sorted_counts(i)/sum(counts));
end
text(0.1, 0.9, stats_text, 'FontSize', 9, 'VerticalAlignment', 'top');

% Distribution comparison
subplot(2, 3, 6);
axis off;
comparison_text = sprintf(['Distribution Analysis:\n\n' ...
                          'Sorted View: Shows Pareto\n' ...
                          'distribution for analysis\n\n' ...
                          'Uniform View: Shows raw\n' ...
                          'sampling distribution\n\n' ...
                          'Both should be similar\n' ...
                          'for uniform sampling']);
text(0.1, 0.9, comparison_text, 'FontSize', 9, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'encounter_types_analysis');
end

function plot_family_statistics(population, families, opts)
% Plot family-wise statistics
figure('Position', [100, 100, 1000, 800]);

% Family counts
subplot(2, 3, 1);
family_counts = zeros(1, length(families));
for i = 1:length(families)
    family_counts(i) = sum(population.baseline_family_idx == i);
end
bar(family_counts);
set(gca, 'XTickLabel', families);
xlabel('Family');
ylabel('Count');
title('Encounters per Family');
grid on;

% Miss distance by family
subplot(2, 3, 2);
boxplot(population.dmiss_km, population.baseline_family_idx);
set(gca, 'XTickLabel', families);
xlabel('Family');
ylabel('Miss Distance [km]');
title('Miss Distance by Family');
grid on;

% Relative velocity by family
subplot(2, 3, 3);
boxplot(population.vrel_kms, population.baseline_family_idx);
set(gca, 'XTickLabel', families);
xlabel('Family');
ylabel('Relative Velocity [km/s]');
title('Relative Velocity by Family');
grid on;

% TCA times by family
subplot(2, 3, 4);
boxplot(population.tca_days, population.baseline_family_idx);
set(gca, 'XTickLabel', families);
xlabel('Family');
ylabel('TCA Time [days]');
title('TCA Times by Family');
grid on;

% Encounter angles by family
subplot(2, 3, 5);
boxplot(population.angles_deg, population.baseline_family_idx);
set(gca, 'XTickLabel', families);
xlabel('Family');
ylabel('Encounter Angle [deg]');
title('Encounter Angles by Family');
grid on;

% Family summary table
subplot(2, 3, 6);
axis off;
summary_text = sprintf('Family Summary:\n\n');
for i = 1:length(families)
    summary_text = sprintf('%s%s: %d encounters\n', ...
                          summary_text, families{i}, family_counts(i));
end
summary_text = sprintf('%s\nTotal: %d encounters', ...
                      summary_text, sum(family_counts));
text(0.1, 0.9, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'family_statistics');
end

function plot_orbital_elements_analysis(states, mu, opts)
% Plot orbital elements analysis for baseline states
figure('Position', [100, 100, 800, 600]);

% Extract position and velocity
r = states(1:3, :);
v = states(4:6, :);

% Calculate orbital elements for each state
n_states = size(states, 2);
sma = zeros(1, n_states);
ecc = zeros(1, n_states);
inc = zeros(1, n_states);

for i = 1:n_states
    % Convert SI->km, call validated car2kep, map outputs back
    rr_km  = r(:, i) / 1000;
    vv_kms = v(:, i) / 1000;
    mu_km3 = mu(i) / 1e9;
    [a_km,e_i,i_i,~,~,~] = car2kep(rr_km, vv_kms, mu_km3);
    sma(i) = a_km * 1000;    % meters
    ecc(i) = e_i;
    inc(i) = i_i;            % radians
end

% Semi-major axis distribution
subplot(2, 3, 1);
histogram(sma/1000, 30, 'Normalization', 'pdf');
xlabel('Semi-major Axis [km]');
ylabel('Probability Density');
title('Semi-major Axis Distribution');
grid on;

% Eccentricity distribution
subplot(2, 3, 2);
histogram(ecc, 30, 'Normalization', 'pdf');
xlabel('Eccentricity');
ylabel('Probability Density');
title('Eccentricity Distribution');
grid on;

% Inclination distribution
subplot(2, 3, 3);
histogram(rad2deg(inc), 30, 'Normalization', 'pdf');
xlabel('Inclination [deg]');
ylabel('Probability Density');
title('Inclination Distribution');
grid on;

% Altitude distribution (approximate): (a - R_moon) in km
% Use configured Moon radius from defaults to avoid hard-coding
cfg_local = LLO_defaults();
R_moon_m = cfg_local.params.R_moon; % meters
altitudes = (sma - R_moon_m) / 1000; % sma in m -> km above lunar surface
subplot(2, 3, 4);
histogram(altitudes, 30, 'Normalization', 'pdf');
xlabel('Approximate Altitude [km]');
ylabel('Probability Density');
title('Altitude Distribution');
grid on;

% Orbital element correlations
subplot(2, 3, 5);
scatter(ecc, altitudes, 10, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Eccentricity');
ylabel('Approximate Altitude [km]');
title('Eccentricity vs Altitude');
grid on;

% Statistics text
subplot(2, 3, 6);
axis off;
stats_text = sprintf(['Orbital Elements:\n\n' ...
                     'SMA: %.1f ± %.1f km\n' ...
                     'Ecc: %.3f ± %.3f\n' ...
                     'Inc: %.1f° ± %.1f°\n' ...
                     'Alt: %.1f ± %.1f km\n\n' ...
                     'Total Orbits: %d'], ...
                     mean(sma)/1000, std(sma)/1000, ...
                     mean(ecc), std(ecc), ...
                     rad2deg(mean(inc)), rad2deg(std(inc)), ...
                     mean(altitudes), std(altitudes), ...
                     n_states);
text(0.1, 0.9, stats_text, 'FontSize', 9, 'VerticalAlignment', 'top');

% Save plot
save_plot(opts, 'orbital_elements_analysis');
end

function stats = generate_population_statistics(population, families)
% Generate comprehensive statistical summary
stats = struct();

% Overall statistics
stats.total_encounters = length(population.dmiss_km);
stats.families = families;
stats.n_families = length(families);

% Miss distance statistics
stats.miss_distance.min = min(population.dmiss_km);
stats.miss_distance.max = max(population.dmiss_km);
stats.miss_distance.mean = mean(population.dmiss_km);
stats.miss_distance.median = median(population.dmiss_km);
stats.miss_distance.std = std(population.dmiss_km);
stats.miss_distance.percentiles = prctile(population.dmiss_km, [2, 25, 50, 75, 98]);

% TCA times statistics
stats.tca_times.min = min(population.tca_days);
stats.tca_times.max = max(population.tca_days);
stats.tca_times.mean = mean(population.tca_days);
stats.tca_times.median = median(population.tca_days);
stats.tca_times.std = std(population.tca_days);
stats.tca_times.percentiles = prctile(population.tca_days, [2, 25, 50, 75, 98]);

% Relative velocity statistics
stats.relative_velocity.min = min(population.vrel_kms);
stats.relative_velocity.max = max(population.vrel_kms);
stats.relative_velocity.mean = mean(population.vrel_kms);
stats.relative_velocity.median = median(population.vrel_kms);
stats.relative_velocity.std = std(population.vrel_kms);
stats.relative_velocity.percentiles = prctile(population.vrel_kms, [2, 25, 50, 75, 98]);

% Encounter angles statistics
stats.encounter_angles.min = min(population.angles_deg);
stats.encounter_angles.max = max(population.angles_deg);
stats.encounter_angles.mean = mean(population.angles_deg);
stats.encounter_angles.median = median(population.angles_deg);
stats.encounter_angles.std = std(population.angles_deg);
stats.encounter_angles.percentiles = prctile(population.angles_deg, [2, 25, 50, 75, 98]);

% Encounter types statistics
[unique_types, ~, ic] = unique(population.encounter_types);
counts = accumarray(ic, 1);
[~, sort_idx] = sort(counts, 'descend');
stats.encounter_types.unique_types = unique_types(sort_idx);
stats.encounter_types.counts = counts(sort_idx);
stats.encounter_types.percentages = 100 * counts(sort_idx) / sum(counts);

% Family-wise statistics
stats.family_stats = struct();
for i = 1:length(families)
    fam_mask = population.baseline_family_idx == i;
    if any(fam_mask)
        stats.family_stats.(families{i}).count = sum(fam_mask);
        stats.family_stats.(families{i}).miss_distance_mean = mean(population.dmiss_km(fam_mask));
        stats.family_stats.(families{i}).vrel_mean = mean(population.vrel_kms(fam_mask));
        stats.family_stats.(families{i}).tca_mean = mean(population.tca_days(fam_mask));
        stats.family_stats.(families{i}).angle_mean = mean(population.angles_deg(fam_mask));
    end
end

% Quality metrics
stats.quality_metrics.generation_success_rate = population.total_generated / population.total_requested;
stats.quality_metrics.average_miss_distance = stats.miss_distance.mean;
stats.quality_metrics.average_relative_velocity = stats.relative_velocity.mean;
end

function print_population_summary(stats)
% Print statistical summary to console
fprintf('\n=== TCA POPULATION STATISTICAL SUMMARY ===\n');
fprintf('Total Encounters: %d\n', stats.total_encounters);
fprintf('Families: %d (%s)\n', stats.n_families, strjoin(stats.families, ', '));

fprintf('\n--- Miss Distance Statistics ---\n');
fprintf('Range: %.2e - %.2e km\n', stats.miss_distance.min, stats.miss_distance.max);
fprintf('Mean: %.2e km, Median: %.2e km, Std: %.2e km\n', ...
        stats.miss_distance.mean, stats.miss_distance.median, stats.miss_distance.std);

fprintf('\n--- TCA Times Statistics ---\n');
fprintf('Range: %.3f - %.3f days\n', stats.tca_times.min, stats.tca_times.max);
fprintf('Mean: %.3f days, Median: %.3f days, Std: %.3f days\n', ...
        stats.tca_times.mean, stats.tca_times.median, stats.tca_times.std);

fprintf('\n--- Relative Velocity Statistics ---\n');
fprintf('Range: %.4f - %.4f km/s\n', stats.relative_velocity.min, stats.relative_velocity.max);
fprintf('Mean: %.4f km/s, Median: %.4f km/s, Std: %.4f km/s\n', ...
        stats.relative_velocity.mean, stats.relative_velocity.median, stats.relative_velocity.std);

fprintf('\n--- Encounter Angles Statistics ---\n');
fprintf('Range: %.1f° - %.1f°\n', stats.encounter_angles.min, stats.encounter_angles.max);
fprintf('Mean: %.1f°, Median: %.1f°, Std: %.1f°\n', ...
        stats.encounter_angles.mean, stats.encounter_angles.median, stats.encounter_angles.std);

fprintf('\n--- Encounter Types (Top 5) ---\n');
for i = 1:min(5, length(stats.encounter_types.unique_types))
    fprintf('%s: %d (%.1f%%)\n', stats.encounter_types.unique_types{i}, ...
            stats.encounter_types.counts(i), stats.encounter_types.percentages(i));
end

fprintf('\n--- Quality Metrics ---\n');
fprintf('Generation Success Rate: %.1f%%\n', 100 * stats.quality_metrics.generation_success_rate);
fprintf('Average Miss Distance: %.2e km\n', stats.quality_metrics.average_miss_distance);
fprintf('Average Relative Velocity: %.4f km/s\n', stats.quality_metrics.average_relative_velocity);

fprintf('\n=== END SUMMARY ===\n\n');
end

function save_plot(opts, plot_name)
% Helper function to save plots
if ~isfield(opts.plot, 'format')
    opts.plot.format = 'png';
end

filename = sprintf('%s_LLO_TCA.%s', plot_name, opts.plot.format);
filepath = fullfile(opts.plot.saveDir, filename);

try
    saveas(gcf, filepath);
    fprintf('Saved plot: %s\n', filepath);
catch ME
    warning('Failed to save plot %s: %s', plot_name, ME.message);
end
end
