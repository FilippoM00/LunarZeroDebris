function [t_days, stats] = sample_tca_times_uniform(min_days, max_days, n, opts)
%SAMPLE_TCA_TIMES_UNIFORM Draw time-to-TCA samples [days] from a uniform distribution.
%
%   [t_days, stats] = sample_tca_times_uniform(min_days, max_days, n)
%   returns n samples uniformly distributed in [min_days, max_days].
%   Optionally plots a histogram with visible edges and saves it.
%
%   Inputs
%     min_days  (scalar): minimum time to TCA in days
%     max_days  (scalar): maximum time to TCA in days (> min_days)
%     n         (integer): number of samples
%     opts      (struct, optional):
%        .seed           RNG seed (default: [])
%        .plot           true/false (default true)
%        .bins           number of histogram bins (default 30)
%        .normalization  'count' (default) | 'pdf' | 'probability' etc.
%        .title          custom title (default auto)
%        .saveTo         file path to save figure (default '')
%
%   Outputs
%     t_days : n-by-1 vector of sampled times to TCA in days
%     stats  : struct with summary statistics (min, max, mean, std, percentiles)
%
%   Example
%     [t, s] = sample_tca_times_uniform(1, 4, 1000, struct('plot',true));

	if nargin < 3
		error('Provide min_days, max_days, and n.');
	end
	if nargin < 4 || isempty(opts), opts = struct(); end

	% Validate
	validateattributes(min_days, {'numeric'},{'scalar','real','finite'} ,mfilename,'min_days',1);
	validateattributes(max_days, {'numeric'},{'scalar','real','finite','>',min_days} ,mfilename,'max_days',2);
	validateattributes(n, {'numeric'},{'scalar','integer','>=',1,'finite'} ,mfilename,'n',3);

	% Defaults
	if ~isfield(opts,'seed'), opts.seed = []; end
	if ~isfield(opts,'plot'), opts.plot = true; end
	if ~isfield(opts,'bins') || isempty(opts.bins), opts.bins = 30; end
	if ~isfield(opts,'normalization') || isempty(opts.normalization), opts.normalization = 'count'; end
	if ~isfield(opts,'title') || isempty(opts.title)
		opts.title = sprintf('Uniform Time-to-TCA Samples (n=%d), [%g,%g] days', n, min_days, max_days);
	end
	if ~isfield(opts,'saveTo'), opts.saveTo = ''; end

	% Seed control
	if ~isempty(opts.seed), rng(opts.seed); end

	% Sample uniformly
	t_days = min_days + (max_days - min_days) * rand(n,1);

	% Stats
	stats = struct();
	stats.min = min(t_days); stats.max = max(t_days);
	stats.mean = mean(t_days); stats.std = std(t_days);
	pr = prctile(t_days, [2 25 50 75 98]);
	stats.p2 = pr(1); stats.p25 = pr(2); stats.p50 = pr(3); stats.p75 = pr(4); stats.p98 = pr(5);

	% Plot
	if opts.plot
		edges = linspace(min_days, max_days, opts.bins+1);
		fig = figure('Name','Time-to-TCA Distribution (Uniform)','Position',[120,120,1000,600]);
		hold on
		histogram(t_days, 'BinEdges', edges, 'Normalization', opts.normalization, ...
				  'FaceColor', [0.3 0.7 0.4], 'FaceAlpha', 0.35, 'EdgeColor', [0.1 0.4 0.2], 'LineWidth', 0.75);
		xline(min_days,'--','min','Color',[0.2 0.2 0.2]);
		xline(max_days,'--','max','Color',[0.2 0.2 0.2]);
		xline(stats.mean,':','mean','Color',[0.2 0.2 0.6]);
		grid on; box on
		xlabel('Time to TCA [days]');
		if strcmpi(opts.normalization,'pdf')
			ylabel('PDF');
		else
			ylabel('Count');
		end
		title(opts.title);
		legend('Samples','Location','best');
		hold off
		if ~isempty(opts.saveTo)
			try, saveas(fig, opts.saveTo); catch, end
		end
	end
end
