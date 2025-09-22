function [types, vrel_kms, angles_deg, meta] = sample_encounter_params_uniform(n, opts)
%SAMPLE_ENCOUNTER_PARAMS_UNIFORM Sample encounter types, relative speeds, and angles.
%
%   [types, vrel_kms, angles_deg, meta] = sample_encounter_params_uniform(n, opts)
%   returns n samples of:
%     - types: cellstr like 'slow_trailing', 'fast_head_on', ...
%     - vrel_kms: relative velocity magnitudes in km/s sampled uniformly within class ranges
%     - angles_deg: approach angles in degrees sampled uniformly within direction ranges
%   Options let you customize speed ranges and angle ranges.
%
%   Inputs
%     n    (integer): number of samples
%     opts (struct, optional):
%       .seed      RNG seed
%       .classes   cellstr of classes to use (default all 15: 5 speeds x 3 dirs)
%       (weights removed; classes sampled uniformly)
%       .speedRanges_kms: struct with fields drift/slow/medium/fast/hyperbolic as [min max] (km/s)
%       .angleRanges_deg: struct with fields trailing/head_on/lateral as [min max] degrees
%       .lateralSymmetric: true to choose around 90 or 270 deg (default true)
%       .plot      true to plot distributions (default false)
%       .saveDir   folder to save plots (default '')
%       .bandWeights struct with fields drift/slow/medium/fast/hyperbolic giving
%                   selection probabilities for speed bands. If provided, they
%                   control how often each band is selected (directions remain
%                   uniform). If omitted and no custom classes list is given,
%                   defaults bias against drift & hyperbolic.
%
%   Outputs
%     types      : n x 1 cellstr of class names
%     vrel_kms   : n x 1 double (km/s)
%     angles_deg : n x 1 double (deg)
%     meta       : struct with the resolved config (classes, ranges)

	if nargin < 1, error('Provide n.'); end
	if nargin < 2 || isempty(opts), opts = struct(); end
	validateattributes(n, {'numeric'},{'scalar','integer','>=',1});

	if ~isfield(opts,'seed'), opts.seed = []; end
	if ~isfield(opts,'plot'), opts.plot = false; end
	if ~isfield(opts,'saveDir'), opts.saveDir = ''; end
	if ~isfield(opts,'lateralSymmetric'), opts.lateralSymmetric = true; end
	if ~isfield(opts,'speedSampling') || isempty(opts.speedSampling), opts.speedSampling = 'uniform'; end
	% Per-band speed sampling defaults: drift log-normal increasing (mode > max),
	% slow/medium/fast log-uniform (flat on log plots), hyperbolic log-normal (decreasing)
	defSSB = struct('drift','lognormal_inc','slow','loguniform','medium','loguniform','fast','loguniform','hyperbolic','lognormal');
	if ~isfield(opts,'speedSamplingByBand') || ~isstruct(opts.speedSamplingByBand)
		ssb = defSSB;
	else
		ssb = defSSB; fn = fieldnames(opts.speedSamplingByBand);
		for i=1:numel(fn), ssb.(fn{i}) = lower(string(opts.speedSamplingByBand.(fn{i}))); end
	end
	% Hyperbolic log-normal parameters
	% Hyperbolic log-normal defaults: place mean kSigma*sigma to the left of ln(min)
	defHyp = struct('sigma_ln', 1.0, 'kSigma', 1.0, 'epsilon', 0.1, 'modeOverride', []);
	if isfield(opts,'hyperbolicLogNormal') && isstruct(opts.hyperbolicLogNormal)
		hypLN = defHyp; fn = fieldnames(opts.hyperbolicLogNormal);
		for i=1:numel(fn), hypLN.(fn{i}) = opts.hyperbolicLogNormal.(fn{i}); end
	else
		hypLN = defHyp;
	end
	% Drift log-normal-increasing parameters (mode to the right of max)
	defDrift = struct('sigma_ln', 1.0, 'kSigma', 1.0, 'modeOverride', []);
	if isfield(opts,'driftLogNormal') && isstruct(opts.driftLogNormal)
		driftLN = defDrift; fn = fieldnames(opts.driftLogNormal);
		for i=1:numel(fn), driftLN.(fn{i}) = opts.driftLogNormal.(fn{i}); end
	else
		driftLN = defDrift;
	end

	% Default ranges (km/s)
	defSpeed.drift      = [0.00005, 0.0005];   % 0.05..0.5 m/s
	defSpeed.slow       = [0.0005, 0.01];      % 0.5..10 m/s
	defSpeed.medium     = [0.01,  0.1];        % 10..100 m/s
	defSpeed.fast       = [0.1,   0.5];        % 0.1..0.5 km/s
	defSpeed.hyperbolic = [0.5,   1.2];        % 0.5..1.2 km/s

	% Default approach angle windows (deg)
	defAngles.trailing = [-45, 45];
	defAngles.head_on  = [135, 225];
	defAngles.lateral  = [45, 135];  % around 90 (or 270)

	speed = defSpeed;
	if isfield(opts,'speedRanges_kms') && isstruct(opts.speedRanges_kms)
		fn = fieldnames(opts.speedRanges_kms);
		for i=1:numel(fn)
			speed.(fn{i}) = opts.speedRanges_kms.(fn{i});
		end
	end
	ang = defAngles;
	if isfield(opts,'angleRanges_deg') && isstruct(opts.angleRanges_deg)
		fn = fieldnames(opts.angleRanges_deg);
		for i=1:numel(fn)
			ang.(fn{i}) = opts.angleRanges_deg.(fn{i});
		end
	end

	% Build class list
	speedNames = {'drift','slow','medium','fast','hyperbolic'};
	dirNames   = {'trailing','head_on','lateral'};
	allClasses = cell(0,1);
	for i=1:numel(speedNames)
		for j=1:numel(dirNames)
			allClasses{end+1,1} = sprintf('%s_%s', speedNames{i}, dirNames{j}); %#ok<AGROW>
		end
	end
	if ~isfield(opts,'classes') || isempty(opts.classes)
		classes = allClasses;
	else
		classes = opts.classes(:);
	end

	% Class selection: uniform across classes

	if ~isempty(opts.seed), rng(opts.seed); end

	% Sample classes
	if isfield(opts,'classes') && ~isempty(opts.classes)
		% Respect a user-provided class list: uniform across those classes
		idx = randi(numel(classes), n, 1);
		types = classes(idx);
	else
		% No class list provided: select speed band with weights, direction uniformly
		bands = {'drift','slow','medium','fast','hyperbolic'};
		% Default band weights (bias against drift & hyperbolic)
		defBW = struct('drift',0.05,'slow',0.30,'medium',0.35,'fast',0.25,'hyperbolic',0.05);
		bw = defBW;
		if isfield(opts,'bandWeights') && isstruct(opts.bandWeights)
			fn = fieldnames(opts.bandWeights);
			for i=1:numel(fn)
				bw.(fn{i}) = opts.bandWeights.(fn{i});
			end
		end
		% Normalize and guard
		wv = [bw.drift, bw.slow, bw.medium, bw.fast, bw.hyperbolic];
		wv(~isfinite(wv) | wv<0) = 0; s = sum(wv);
		if s <= 0, wv = [0.10, 0.30, 0.35, 0.20, 0.05]; s = sum(wv); end
		wv = wv / s;
		% Cumulative for sampling
		cumw = cumsum(wv);
		r = rand(n,1);
		bandIdx = arrayfun(@(x) find(x <= cumw, 1, 'first'), r);
		dirNames = {'trailing','head_on','lateral'};
		dirIdx = randi(numel(dirNames), n, 1);
		types = cell(n,1);
		for k=1:n
			types{k} = sprintf('%s_%s', bands{bandIdx(k)}, dirNames{dirIdx(k)});
		end
	end

	% Sample per class
	vrel_kms = zeros(n,1);
	angles_deg = zeros(n,1);
	for k=1:n
		t = string(types{k});
		parts = strsplit(t, '_');
		if numel(parts) < 2
			sp = 'medium'; dr = 'lateral';
		else
			sp = parts{1};
			% Rejoin remaining pieces to support directions like 'head_on'
			dr = strjoin(parts(2:end), '_');
		end
		% Fallback safeguards if keys are missing
		if ~isfield(speed, sp), sp = 'medium'; end
		if ~isfield(ang, dr), dr = 'lateral'; end
		% speed (per-band sampling)
		rngV = speed.(sp);
		% Determine mode for this band: per-band overrides global; fallback to global
		if isfield(ssb, sp)
			modeSp = char(ssb.(sp));
		else
			modeSp = opts.speedSampling; % fallback
		end
		switch lower(modeSp)
			case {'loguniform','log-uniform'}
				vrel_kms(k) = 10.^(log10(rngV(1)) + (log10(rngV(2)) - log10(rngV(1))) * rand());
			case {'lognormal_inc','lognormal-increasing'}
				% Draw from truncated log-normal with mode beyond max (increasing across [min,max])
				vrel_kms(k) = draw_lognormal_truncated_increasing(rngV(1), rngV(2), driftLN);
			case 'lognormal'
				% Draw from truncated log-normal decreasing over [min,max]
				vrel_kms(k) = draw_lognormal_truncated(rngV(1), rngV(2), hypLN);
			otherwise
				% linear uniform
				vrel_kms(k) = rngV(1) + (rngV(2) - rngV(1)) * rand();
		end
		% angle (uniform within direction windows)
		if strcmp(dr,'lateral') && opts.lateralSymmetric
			center = 90; if rand() < 0.5, center = 270; end
			spread = ang.lateral - 90; % e.g., [45,135]-90 = [-45,45]
			a = center + (spread(1) + (spread(2)-spread(1))*rand());
			angles_deg(k) = a;
		else
			rngA = ang.(dr);
			angles_deg(k) = rngA(1) + (rngA(2) - rngA(1)) * rand();
		end
		% wrap to [0,360)
		angles_deg(k) = mod(angles_deg(k), 360);
	end

	% Meta
	meta = struct();
	meta.classes = classes;
	meta.speedRanges_kms = speed; meta.angleRanges_deg = ang;
	meta.speedSampling = opts.speedSampling;
	meta.speedSamplingByBand = ssb;
	% Report effective log-normal params based on current ranges
	try
		rngH = speed.hyperbolic; % [min,max]
		[mu_eff, mode_eff] = local_hyp_lognormal_params(rngH(1), hypLN);
		meta.hyperbolicLogNormal = struct('sigma_ln', hypLN.sigma_ln, 'kSigma', hypLN.kSigma, 'epsilon', hypLN.epsilon, 'modeOverride', hypLN.modeOverride, 'mu_ln', mu_eff, 'mode_ln', mode_eff);
	catch
		meta.hyperbolicLogNormal = hypLN;
	end
	try
		rngD = speed.drift; % [min,max]
		[mu_d, mode_d] = local_drift_lognormal_params(rngD(2), driftLN);
		meta.driftLogNormal = struct('sigma_ln', driftLN.sigma_ln, 'kSigma', driftLN.kSigma, 'modeOverride', driftLN.modeOverride, 'mu_ln', mu_d, 'mode_ln', mode_d);
	catch
		meta.driftLogNormal = driftLN;
	end
	if ~(isfield(opts,'classes') && ~isempty(opts.classes))
		meta.bandWeights = struct('drift',bw.drift,'slow',bw.slow,'medium',bw.medium,'fast',bw.fast,'hyperbolic',bw.hyperbolic);
	end

	% Optional plots
	if opts.plot
		if ~isempty(opts.saveDir) && ~exist(opts.saveDir,'dir'), mkdir(opts.saveDir); end
		% 1) Counts per class
		[uC,~,ic] = unique(types);
		cnt = accumarray(ic, 1);
		fig1 = figure('Name','Encounter Types','Position',[80,80,1000,500]);
		bar(categorical(uC), cnt, 'FaceColor',[0.4 0.6 0.9]); grid on; box on
		title('Counts by Encounter Type'); ylabel('Count');
		if ~isempty(opts.saveDir)
			try, saveas(fig1, fullfile(opts.saveDir,'encounter_types_counts.png')); catch, end
		end
		% 2) v_rel histogram
		fig2 = figure('Name','Relative Speeds','Position',[100,120,1000,500]);
		% Build per-band log-spaced edges to avoid bins straddling boundaries
		bands = {'drift','slow','medium','fast','hyperbolic'};
		binsPerBand = 10;
		edgesV = [];
		for bi = 1:numel(bands)
			br = speed.(bands{bi});
			bmin = max(br(1), 1e-8);
			bmax = max(br(2), bmin*(1+1e-6));
			e = logspace(log10(bmin), log10(bmax), binsPerBand+1);
			if ~isempty(edgesV) && abs(e(1) - edgesV(end)) < 1e-16
				e = e(2:end); % avoid duplicate edge at junction
			end
			edgesV = [edgesV, e]; %#ok<AGROW>
		end
		histogram(vrel_kms, 'BinEdges', edgesV, 'FaceColor',[0.3 0.7 0.4], 'FaceAlpha',0.4, 'EdgeColor',[0.1 0.4 0.2], 'LineWidth',0.9);
		set(gca,'XScale','log'); grid on; box on
		xlabel('Relative speed [km/s]'); ylabel('Count'); title('v_{rel} distribution (log scale)');
		hold on
		% Band boundary markers for clarity
		for bi = 1:numel(bands)
			br = speed.(bands{bi});
			xline(br(1), ':', 'Color',[0.2 0.2 0.2], 'LineWidth',0.8);
			xline(br(2), ':', 'Color',[0.2 0.2 0.2], 'LineWidth',0.8);
		end
		hold off
		if ~isempty(opts.saveDir)
			try, saveas(fig2, fullfile(opts.saveDir,'vrel_hist.png')); catch, end
		end
		% 3) Angle histogram
		fig3 = figure('Name','Approach Angles','Position',[120,160,1000,500]);
		edgesA = 0:5:360;
		histogram(angles_deg, 'BinEdges', edgesA, 'FaceColor',[0.9 0.6 0.3], 'FaceAlpha',0.35, 'EdgeColor',[0.6 0.3 0.1], 'LineWidth',0.75);
		grid on; box on; xlim([0,360]);
		xlabel('Approach angle [deg]'); ylabel('Count'); title('Approach angle distribution');
		if ~isempty(opts.saveDir)
			try, saveas(fig3, fullfile(opts.saveDir,'angle_hist.png')); catch, end
		end
	end
end

function v = draw_lognormal_truncated(minV, maxV, hypLN)
%DRAW_LOGNORMAL_TRUNCATED Sample from log-normal truncated to [minV,maxV].
% The parameters are configured so the mode is below minV, yielding a
% decreasing density on [minV, +inf).
	sigma = hypLN.sigma_ln;
	if ~isfinite(minV) || ~isfinite(maxV) || minV<=0 || maxV<=0 || maxV<=minV
		% Guard: fallback to linear uniform
		v = minV + (maxV - minV) * rand();
		return
	end
	if isfield(hypLN,'modeOverride') && ~isempty(hypLN.modeOverride) && isfinite(hypLN.modeOverride) && hypLN.modeOverride>0
		mu = log(hypLN.modeOverride) + sigma^2; % mode = exp(mu - sigma^2)
	else
		% Place log-space mean to the left of ln(minV) to ensure decreasing trend
		kSig = hypLN.kSigma; if ~isfinite(kSig) || kSig < 0, kSig = 1.0; end
		mu = log(minV) - kSig * sigma;
		% Keep legacy epsilon: also pushes the mode further below min if desired
	% epsilon parameter is retained for backward compatibility; not used here
	end
	% Rejection sampling
	maxIters = 200;
	for it=1:maxIters
		z = mu + sigma*randn(); v = exp(z);
		if v>=minV && v<=maxV
			return
		end
	end
	% Fallback clamp if rejection failed
	v = min(max(v, minV), maxV);
end

function [mu_ln, mode_val] = local_hyp_lognormal_params(minV, hypLN)
%LOCAL_HYP_LOGNORMAL_PARAMS Compute effective mu and mode for reporting
	sigma = hypLN.sigma_ln;
	if isfield(hypLN,'modeOverride') && ~isempty(hypLN.modeOverride) && isfinite(hypLN.modeOverride) && hypLN.modeOverride>0
		mu_ln = log(hypLN.modeOverride) + sigma^2;
		mode_val = hypLN.modeOverride;
	else
		kSig = hypLN.kSigma; if ~isfinite(kSig) || kSig < 0, kSig = 1.0; end
		mu_ln = log(minV) - kSig * sigma;
		mode_val = exp(mu_ln - sigma^2);
	end
end

function v = draw_lognormal_truncated_increasing(minV, maxV, driftLN)
%DRAW_LOGNORMAL_TRUNCATED_INCREASING Log-normal with mode beyond maxV; truncated to [minV,maxV]
	sigma = driftLN.sigma_ln;
	if ~isfinite(minV) || ~isfinite(maxV) || minV<=0 || maxV<=0 || maxV<=minV
		v = minV + (maxV - minV) * rand();
		return
	end
	if isfield(driftLN,'modeOverride') && ~isempty(driftLN.modeOverride) && isfinite(driftLN.modeOverride) && driftLN.modeOverride>0
		mu = log(driftLN.modeOverride) + sigma^2; % mode = exp(mu - sigma^2)
	else
		kSig = driftLN.kSigma; if ~isfinite(kSig) || kSig < 0, kSig = 1.0; end
		% Place mean to the right of ln(maxV)
		mu = log(maxV) + kSig * sigma;
	end
	% Rejection sampling
	maxIters = 200;
	for it=1:maxIters
		z = mu + sigma*randn(); v = exp(z);
		if v>=minV && v<=maxV
			return
		end
	end
	v = min(max(v, minV), maxV);
end

function [mu_ln, mode_val] = local_drift_lognormal_params(maxV, driftLN)
%LOCAL_DRIFT_LOGNORMAL_PARAMS Effective mu and mode for drift (increasing)
	sigma = driftLN.sigma_ln;
	if isfield(driftLN,'modeOverride') && ~isempty(driftLN.modeOverride) && isfinite(driftLN.modeOverride) && driftLN.modeOverride>0
		mu_ln = log(driftLN.modeOverride) + sigma^2;
		mode_val = driftLN.modeOverride;
	else
		kSig = driftLN.kSigma; if ~isfinite(kSig) || kSig < 0, kSig = 1.0; end
		mu_ln = log(maxV) + kSig * sigma;
		mode_val = exp(mu_ln - sigma^2);
	end
end
