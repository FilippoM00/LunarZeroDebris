function [samples_km, params] = sample_miss_distances_lognormal(min_km, max_km, n, opts)
%SAMPLE_MISS_DISTANCES_LOGNORMAL Draw miss distances [km] from a log-normal distribution.
%
%   [samples_km, params] = sample_miss_distances_lognormal(min_km, max_km, n)
%   fits a log-normal distribution such that min_km and max_km correspond to
%   the lower/upper coverage quantiles (default 2% and 98%), samples n values,
%   and plots the histogram with the theoretical PDF overlay.
%
%   Inputs
%     min_km (scalar > 0): minimum plausible miss distance in km
%     max_km (scalar > min_km): maximum plausible miss distance in km
%     n      (positive integer): number of samples
%     opts   (struct, optional):
%        .coverage       [qL,qH] quantiles that map to [min_km,max_km] (default [0.02 0.98])
%        .seed           RNG seed (default: [])
%        .mode           'clamp' (default) | 'resample' | 'truncate' (alias of 'resample')
%        .plot           true/false (default true)
%        .bins           number of histogram bins (default 30)
%        .logX           true to use logarithmic X axis (default false)
%        .normalization  'pdf' (default) | 'count' | 'probability' etc. (histogram Normalization)
%        .title          custom plot title (default auto)
%        .saveTo         file path to save figure (e.g., 'miss_dist_hist.png') (default '')
%
%   Outputs
%     samples_km : n-by-1 sampled miss distances in km
%     params     : struct with fitted log-normal parameters and stats
%                  .mu, .sigma (log-space parameters)
%                  .coverage, .zQuantiles, .median_km, .mean_km, .std_km
%                  .mode, .bins, .logX, .normalization, .seedUsed
%
%   Notes
%     - The log-normal is fitted in log-space so that the chosen coverage
%       quantiles match the provided [min_km, max_km]. By default, ~96% of
%       the probability mass lies within [min_km, max_km].
%     - 'clamp' mode softly bounds samples to [min_km, max_km]; 'resample'
%       repeats sampling until all values fall inside the range.
%
%   Example
%     [x, p] = sample_miss_distances_lognormal(0.001, 3.0, 2000, struct('logX',true));
%     % x are samples in km; p.mu/p.sigma define the log-normal
%
%   This function is self-contained and has no toolbox dependencies.

	if nargin < 3
		error('Provide min_km, max_km, and n.');
	end
	if nargin < 4 || isempty(opts)
		opts = struct();
	end

	% Validate inputs
	validateattributes(min_km, {'numeric'},{'scalar','real','positive','finite'} ,mfilename,'min_km',1);
	validateattributes(max_km, {'numeric'},{'scalar','real','finite','>',min_km} ,mfilename,'max_km',2);
	validateattributes(n, {'numeric'},{'scalar','integer','>=',1,'finite'} ,mfilename,'n',3);

	% Defaults
	if ~isfield(opts,'coverage') || ~isnumeric(opts.coverage) || numel(opts.coverage)~=2
		opts.coverage = [0.02 0.98];
	end
	qL = max(min(opts.coverage(1), 0.999999), 1e-6);
	qH = max(min(opts.coverage(2), 0.999999), 1e-6);
	if qL >= qH
		qL = 0.02; qH = 0.98;
	end
	if ~isfield(opts,'seed'), opts.seed = []; end
	if ~isfield(opts,'mode') || isempty(opts.mode), opts.mode = 'clamp'; end
	if ~isfield(opts,'plot'), opts.plot = true; end
	if ~isfield(opts,'bins') || isempty(opts.bins), opts.bins = 30; end
	if ~isfield(opts,'logX'), opts.logX = false; end
	if ~isfield(opts,'normalization') || isempty(opts.normalization), opts.normalization = 'pdf'; end
	if ~isfield(opts,'title') || isempty(opts.title)
		opts.title = sprintf('Log-normal Miss Distance Samples (n=%d), [%g,%g] km', n, min_km, max_km);
	end
	if ~isfield(opts,'saveTo'), opts.saveTo = ''; end

	% Seed control
	if ~isempty(opts.seed)
		rng(opts.seed);
	end

	% Solve for log-normal parameters given quantiles
	z = @(p) -sqrt(2) * erfcinv(2*p); % norminv using base functions
	zL = z(qL); zH = z(qH);
	ln_min = log(min_km); ln_max = log(max_km);
	sigma = (ln_max - ln_min) / (zH - zL);
	mu    = ln_min - sigma * zL;

	% Sample
	raw = exp(mu + sigma*randn(n,1));
	switch lower(string(opts.mode))
		case 'clamp'
			samples_km = min(max(raw, min_km), max_km);
		case {'resample','truncate'}
			samples_km = zeros(n,1);
			i = 1;
			max_iters = 10000; it = 0;
			while i <= n && it < max_iters
				it = it + 1;
				x = exp(mu + sigma*randn(n - i + 1, 1));
				in = (x >= min_km) & (x <= max_km);
				k = sum(in);
				if k > 0
					samples_km(i:i+k-1) = x(in);
					i = i + k;
				end
			end
			if i <= n
				% Fallback to clamping remaining to guarantee size
				r = exp(mu + sigma*randn(n - i + 1, 1));
				samples_km(i:end) = min(max(r, min_km), max_km);
			end
		otherwise
			warning('Unknown mode ''%s''; using clamp.', string(opts.mode));
			samples_km = min(max(raw, min_km), max_km);
	end

	% Stats
	median_km = exp(mu);
	mean_km   = exp(mu + 0.5*sigma^2);
	std_km    = sqrt((exp(sigma^2)-1) * exp(2*mu + sigma^2));

	params = struct();
	params.mu = mu; params.sigma = sigma;
	params.coverage = [qL, qH]; params.zQuantiles = [zL, zH];
	params.median_km = median_km; params.mean_km = mean_km; params.std_km = std_km;
	params.mode = char(opts.mode); params.bins = opts.bins; params.logX = opts.logX; params.normalization = opts.normalization;
	params.seedUsed = opts.seed;

	% Plot
	if opts.plot
		if opts.logX
			edges = logspace(log10(min_km), log10(max_km), opts.bins+1);
			xPlot = logspace(log10(min_km), log10(max_km), 400);
		else
			edges = linspace(min_km, max_km, opts.bins+1);
			xPlot = linspace(min_km, max_km, 400);
		end
		fig = figure('Name','Miss Distance Distribution','Position',[120,120,1000,600]);
		hold on
	histogram(samples_km, 'BinEdges', edges, 'Normalization', opts.normalization, ...
		  'FaceColor', [0.2 0.6 0.9], 'FaceAlpha', 0.35, 'EdgeColor', [0.1 0.3 0.5], 'LineWidth', 0.75);
		% Theoretical PDF
		pdf_vals = (1./(xPlot * sigma * sqrt(2*pi))) .* exp(-0.5*((log(xPlot)-mu)/sigma).^2);
		if strcmpi(opts.normalization,'pdf')
			plot(xPlot, pdf_vals, 'k-', 'LineWidth', 1.5, 'DisplayName','Log-normal PDF');
		else
			% Convert pdf to counts approximately by multiplying by total area under histogram
			% Using average bin width as scale (rough visual overlay)
			bw = mean(diff(edges));
			plot(xPlot, pdf_vals * numel(samples_km) * bw, 'k-', 'LineWidth', 1.5, 'DisplayName','Log-normal PDF (scaled)');
		end
		% Markers
		xline(min_km,'--','q_{low}','Color',[0.2 0.2 0.2]);
		xline(max_km,'--','q_{high}','Color',[0.2 0.2 0.2]);
		xline(median_km,':','median','Color',[0.1 0.5 0.1]);
		if opts.logX
			set(gca,'XScale','log');
		end
		grid on; box on
		xlabel('Miss distance [km]');
		if strcmpi(opts.normalization,'pdf')
			ylabel('PDF');
		else
			ylabel('Count');
		end
		title(opts.title);
		legend('Location','best');
		hold off
		if ~isempty(opts.saveTo)
			try
				saveas(fig, opts.saveTo);
			catch
				% ignore save errors
			end
		end
	end
end
