function [population, meta] = create_tca_population(total, opts)
%CREATE_TCA_POPULATION Entry point to generate a TCA population by count.
%
%   [population, meta] = create_tca_population(total, opts)
%
%   Inputs
%     total (integer > 0): total number of close approaches to generate
%     opts  (struct, optional): configuration options (to be defined)
%
%   Outputs
%     population : struct (to be defined)
%     meta       : struct (to be defined)
%
%   Notes
%     - Implementation will be developed from scratch in subsequent steps.
%     - Secondary states at epoch are constructed to satisfy TCA conditions:
%       r_rel · v_rel = 0 and |r_rel| = d_miss. An optional isotropic 3D
%       randomization rotates r_rel around v_rel to ensure uniform coverage
%       in the perpendicular plane without altering v_rel.

% ======================================================================
% Step 1: load, sample, and plot a fixed pool of primary orbits from
% Database_Orbits using verified defaults. The pool size is independent of
% the requested number of cases: by default the sampler selects 100 orbits
% per family with its built-in policies (no overrides here).
% ======================================================================

	% --- Validate inputs ---
	if nargin < 1
		error('create_tca_population:arg', 'Provide total number of close approaches (e.g., 100).');
	end
	validateattributes(total, {'numeric'},{'scalar','integer','>=',1,'finite'}, mfilename,'total',1);
	if nargin < 2 || isempty(opts)
		opts = struct();
	end

	% Initialize outputs to avoid unset warnings in early stages
	population = struct();
	meta = struct();

	% --- Resolve repository and database paths ---
	here = fileparts(mfilename('fullpath'));
	repo = fileparts(fileparts(here));  % repo root
	if isfield(opts, 'dbRoot') && ~isempty(opts.dbRoot)
		dbRoot = opts.dbRoot;
	else
		dbRoot = fullfile(repo, 'Database_Orbits');
	end
	if exist(dbRoot, 'dir') ~= 7
		error('create_tca_population:dbRoot', 'Database_Orbits not found at: %s', dbRoot);
	end

	% Ensure repository code (CR3BP, utilities) and verified codes are on path
	addpath(genpath(repo));

	% --- Sample primaries with verified defaults (100 per family) ---
	% Use same defaults as validation smoke (seed 12345 unless overridden)
	primariesOpts = struct();
	% Primaries: default seed 12345 (as in sampler); allow global opts.seed override
	if isfield(opts, 'seed') && ~isempty(opts.seed)
		primariesOpts.seed = opts.seed;
	end
	primariesSample = sample_primaries_from_db(dbRoot, primariesOpts);
	% Families list for convenience
	families = {'DRO','L1_halo','L2_halo','L1_lyapunov','L2_lyapunov'};

	% --- Package outputs ---
	population.primaries = primariesSample;        % sampled primaries by family
	population.total_requested = total;
	population.families = families;

	meta.stage = 'primaries_sampling';
	meta.dbRoot = dbRoot;
	if isfield(primariesSample, 'counts_requested')
		meta.counts_requested = primariesSample.counts_requested;
	end
	if isfield(primariesSample, 'counts_actual')
		meta.counts_actual = primariesSample.counts_actual;
	end
	meta.timestamp = datestr(now);
	meta.notes = 'Sampled primaries using verified defaults (100 per family).';


	% ======================================================================
	% Step 2: Using the requested input 'total', sample miss distances [km]
	%         and time-to-TCA [days] with verified samplers. These will
	%         define the encounter population size for later steps.
	% 		  Prepare per-family stratified sampling containers
	% ======================================================================

	% Defaults for miss distance and time bounds
	if ~isfield(opts,'miss') || ~isstruct(opts.miss), opts.miss = struct(); end
	% Match smoke test default lower bound (1 meter = 1e-6 km)
	if ~isfield(opts.miss,'min_km') || isempty(opts.miss.min_km), opts.miss.min_km = 1e-6; end
	if ~isfield(opts.miss,'max_km') || isempty(opts.miss.max_km), opts.miss.max_km = 3.0; end
	if ~isfield(opts,'time') || ~isstruct(opts.time), opts.time = struct(); end
	if ~isfield(opts.time,'min_days') || isempty(opts.time.min_days), opts.time.min_days = 1.0; end
	if ~isfield(opts.time,'max_days') || isempty(opts.time.max_days), opts.time.max_days = 4.0; end

	% Units and isotropy options
	if ~isfield(opts,'units') || ~isstruct(opts.units), opts.units = struct(); end
	if ~isfield(opts.units,'LU_km') || isempty(opts.units.LU_km), opts.units.LU_km = 384400; end % Earth-Moon distance
	if ~isfield(opts.units,'TU_days') || isempty(opts.units.TU_days)
		% TU = T_sidereal / (2*pi). Default to Moon sidereal period.
		opts.units.TU_days = 27.321661 / (2*pi);
	end
	if ~isfield(opts,'isotropic') || ~isstruct(opts.isotropic), opts.isotropic = struct(); end
	if ~isfield(opts.isotropic,'enable') || isempty(opts.isotropic.enable), opts.isotropic.enable = true; end
	if ~isfield(opts.isotropic,'seed'), opts.isotropic.seed = []; end
	if ~isfield(opts.isotropic,'mode') || isempty(opts.isotropic.mode), opts.isotropic.mode = 'around_vrel'; end

	population.n_encounters = total;
	enc_types_all = cell(total,1);
	vrel_all = zeros(total,1);
	angles_all = zeros(total,1);
	dmiss_all = zeros(total,1);
	tca_all = zeros(total,1);
	meta.encounter = struct('by_family', struct());

	% ======================================================================
	% Step 4: Define number of encounters per family (no assignment yet)
	% ======================================================================
	n_families = numel(families);
	n_per_family = floor(total / n_families);
	remainder = rem(total, n_families);
	counts_per_family = repmat(n_per_family, n_families, 1);
	if remainder > 0
		counts_per_family(1:remainder) = counts_per_family(1:remainder) + 1;
	end
	meta.encounter_counts = struct('counts_per_family', counts_per_family, 'families', {families});

	% ======================================================================
	% Step 5: For each family, sample a primary IC, propagate a random time,
	%         and stratify encounter params, TCA times, and miss distances
	% ======================================================================
	ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
	baseline_states = zeros(6, total);
	baseline_mu = zeros(1, total);
	baseline_family_idx = zeros(1, total);
	baseline_primary_idx = zeros(1, total);
	baseline_tau = zeros(1, total); % normalized TU
	k = 0;
	for fi = 1:n_families
		fam = families{fi};
		src = population.primaries.(fam);
		m = counts_per_family(fi);
		if src.count == 0 && m > 0
			error('create_tca_population:noOrbitsInFamily', 'Family %s has no sampled primaries.', fam);
		end
		% Per-family encounter parameters
		encOpts = struct();
		if isfield(opts,'seed') && ~isempty(opts.seed), encOpts.seed = opts.seed; end
		if isfield(opts,'encounter') && isstruct(opts.encounter)
			fn = fieldnames(opts.encounter);
			for ii=1:numel(fn), encOpts.(fn{ii}) = opts.encounter.(fn{ii}); end
		end
		[types_f, vrel_f, angles_f, enc_meta_f] = sample_encounter_params_uniform(m, encOpts);
		meta.encounter.by_family.(fam) = enc_meta_f;

		% Per-family TCA times
		timeOpts = struct('seed',123,'plot',false);
		[tca_f, time_meta_f] = sample_tca_times_uniform(opts.time.min_days, opts.time.max_days, m, timeOpts);
		if ~isfield(meta,'time_by_family'), meta.time_by_family = struct(); end
		meta.time_by_family.(fam) = time_meta_f;

		% Per-family miss distances
		[dmiss_f, miss_meta_f] = sample_miss_distances_lognormal(opts.miss.min_km, opts.miss.max_km, m, struct('seed',42,'plot',false,'coverage',[0.02 0.98],'mode','resample'));
		if ~isfield(meta,'miss_by_family'), meta.miss_by_family = struct(); end
		meta.miss_by_family.(fam) = miss_meta_f;

		for j = 1:m
			idx = randi(src.count); % pick one of the ~100 available
			x0 = src.states(:, idx);
			mu_i = src.mu(idx);
			T = src.T(idx);
			tau = rand() * T; % random time in [0, T)
			try
				[~, X] = ode113(@(t,x) CR3BPunpert(t,x,mu_i), [0, tau], x0(:)', ode_options);
			catch
				[~, X] = ode113(@(t,x) CR3BPunpert(t,x,mu_i), [0, tau], x0(:)', ode_options);
			end
			x_tau = X(end, :)';
			k = k + 1;
			baseline_states(:, k) = x_tau;
			baseline_mu(k) = mu_i;
			baseline_family_idx(k) = fi;
			baseline_primary_idx(k) = idx;
			baseline_tau(k) = tau;
			% Attach stratified encounter/time/miss aligned to k
			enc_types_all{k} = types_f{j};
			vrel_all(k) = vrel_f(j);
			angles_all(k) = angles_f(j);
			dmiss_all(k) = dmiss_f(j);
			tca_all(k) = tca_f(j);
		end
	end
	% Attach baseline outputs
	population.baseline_states = baseline_states;
	population.baseline_mu = baseline_mu;
	population.baseline_family_idx = baseline_family_idx;
	population.baseline_primary_idx = baseline_primary_idx;
	population.baseline_tau = baseline_tau;
	% Attach stratified vectors
	population.encounter_types = enc_types_all;
	population.vrel_kms = vrel_all;
	population.angles_deg = angles_all;
	population.dmiss_km = dmiss_all;
	population.tca_days = tca_all;

	% ======================================================================
	% Step 6: Build secondary states at epoch satisfying TCA constraints.
	%         v_rel has given magnitude and in-plane direction; r_rel is
	%         orthogonal to v_rel with magnitude equal to d_miss. Optionally
	%         randomize r_rel azimuth in the plane perpendicular to v_rel.
	% ======================================================================

	% Prepare RNG for isotropic azimuths
	if ~isempty(opts.isotropic.seed)
		rng(opts.isotropic.seed);
	end

	N = total;
	secondary_states = zeros(6, N);
	r_rel_all = zeros(3, N);
	v_rel_all = zeros(3, N);

	LU_km = opts.units.LU_km;
	TU_sec = opts.units.TU_days * 86400.0;

	for k2 = 1:N
		xp = baseline_states(:, k2);
		rp = xp(1:3);
		vp = xp(4:6);

		% Local encounter frame
		vt = norm(vp);
		if vt < 1e-12
			% Fallback: build a frame from rp
			e_r = rp / max(norm(rp), eps);
			e_t = [0;1;0];
			if abs(dot(e_t,e_r)) > 0.99
				e_t = [1;0;0];
			end
			e_t = e_t - dot(e_t,e_r)*e_r; e_t = e_t / max(norm(e_t), eps);
			e_n = cross(e_t, e_r); e_n = e_n / max(norm(e_n), eps);
			e_l = cross(e_n, e_t); e_l = e_l / max(norm(e_l), eps);
		else
			e_t = vp / vt;                            % tangential
			e_r = rp / max(norm(rp), eps);            % radial
			e_n = cross(e_t, e_r); e_n = e_n / max(norm(e_n), eps); % normal
			e_l = cross(e_n, e_t); e_l = e_l / max(norm(e_l), eps); % lateral (in-plane)
		end

		% v_rel direction from in-plane angle
		theta = angles_all(k2); % degrees
		u_v = cosd(theta)*e_t + sind(theta)*e_l;
		u_v = u_v / max(norm(u_v), eps);

		% Magnitudes in LU/TU and LU
		vrel_mag = vrel_all(k2) * (TU_sec / LU_km); % (km/s) * (s/TU) / (km/LU)
		dmiss_mag = dmiss_all(k2) / LU_km;          % km / (km/LU) = LU

		v_rel_vec = vrel_mag * u_v;

		% Build orthonormal basis of plane perpendicular to u_v
		e_b1 = cross(u_v, e_n); if norm(e_b1) < 1e-12, e_b1 = cross(u_v, e_t); end
		e_b1 = e_b1 / max(norm(e_b1), eps);
		e_b2 = cross(u_v, e_b1); e_b2 = e_b2 / max(norm(e_b2), eps);

		if opts.isotropic.enable
			phi = 2*pi*rand();
			b_hat = cos(phi)*e_b1 + sin(phi)*e_b2;
		else
			% Deterministic default orientation (project radial onto perp plane)
			b_hat = e_b1;
		end
		r_rel_vec = dmiss_mag * b_hat;

		% Compose secondary state at epoch
		xs = [rp + r_rel_vec; vp + v_rel_vec];
		secondary_states(:, k2) = xs;
		r_rel_all(:, k2) = r_rel_vec;
		v_rel_all(:, k2) = v_rel_vec;
	end

	% Attach to outputs
	population.secondary_states = secondary_states;
	population.r_rel = r_rel_all;
	population.v_rel = v_rel_all;

	% Simple quality checks (orthogonality and magnitude)
	try
		dd = zeros(1,N); mm = zeros(1,N);
		for q=1:N
			dd(q) = abs(dot(r_rel_all(:,q), v_rel_all(:,q)));
			mm(q) = abs(norm(r_rel_all(:,q))*LU_km - dmiss_all(q)); % km error
		end
		meta.validation = struct('orthogonality_max', max(dd), 'miss_distance_km_max_error', max(mm));
	catch
	end

end

