% VISUALIZE_POPULATION_DISTRIBUTIONS
% Generate LLO primaries using defaults (no propagation) and visualize
% the distributions used by PoC analysis: altitude, eccentricity, inclination.

clear; clc; close all;

% Ensure Population_generation root on path
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
popRoot = fileparts(thisDir); % ..\Population_generation
addpath(genpath(popRoot));

cfg = LLO_defaults();
params = cfg.params; %#ok<NASGU>

% Families to generate and sample size per family
if isfield(cfg, 'families_to_generate')
    families = cfg.families_to_generate;
else
    families = {'Circular','Eccentric','Polar','Frozen','General','HighEccentric'};
end
if isfield(cfg, 'n_primaries_per_family')
    N = cfg.n_primaries_per_family;
else
    N = 1000; % larger sample for smoother histograms
end

fprintf('Generating %d samples per family (%d total) using LLO_defaults...\n', N, N*numel(families));
all_primaries = [];
for i = 1:numel(families)
    fam = families{i};
    prim = generate_LLO_primaries_MCI(N, fam, []); % use defaults internally
    if isempty(all_primaries)
        all_primaries = prim;
    else
        all_primaries = [all_primaries; prim]; %#ok<AGROW>
    end
end

fprintf('Total rows: %d\n', height(all_primaries));

% Prepare output folder
outDir = fullfile(thisDir, 'Media_MCI_Analysis');
if ~exist(outDir, 'dir'); mkdir(outDir); end
ts = datestr(now,'yyyy-mm-dd_HHMMSS');

% Variables to plot
vars = { 'Altitude_km', 'Eccentricity', 'Inclination_deg' };
labels = { 'Altitude (km)', 'Eccentricity (-)', 'Inclination (deg)' };

% Plot distributions by family
fam_list = unique(all_primaries.Family);
colors = lines(numel(fam_list));

for vi = 1:numel(vars)
    var = vars{vi};
    figure('Name', ['Distribution: ' var], 'NumberTitle', 'off', 'Position', [100,100,1200,700]);
    tiledlayout(2, ceil(numel(fam_list)/2), 'Padding','compact','TileSpacing','compact');
    for fi = 1:numel(fam_list)
        nexttile; hold on; grid on;
        fam = fam_list{fi};
        sel = strcmp(all_primaries.Family, fam);
        data = all_primaries.(var)(sel);
        % Choose bins adaptively
        nbins = max(20, round(sqrt(sum(sel))));
        histogram(data, nbins, 'FaceColor', colors(fi,:), 'EdgeColor', 'none');
        xlabel(labels{vi}); ylabel('Count');
        title(sprintf('%s (N=%d)', fam, sum(sel)));

        % Overlay default ranges where applicable
        if strcmp(var, 'Altitude_km')
            if strcmpi(fam,'Circular') || strcmpi(fam,'Polar')
                rng = cfg.alt_range;
            elseif strcmpi(fam,'Eccentric')
                rng = cfg.ecc_peri_alt_range;
            elseif startsWith(fam,'Frozen')
                ff = fam;
                if isfield(cfg.frozen_families, ff) && isfield(cfg.frozen_families.(ff),'alt_range_km')
                    rng = cfg.frozen_families.(ff).alt_range_km;
                else
                    rng = [];
                end
            else
                rng = [];
            end
            if ~isempty(rng)
                yl = ylim;
                plot([rng(1) rng(1)], yl, 'k--', 'LineWidth', 1);
                plot([rng(2) rng(2)], yl, 'k--', 'LineWidth', 1);
                ylim(yl);
                legend({'Samples','Default range'}, 'Location','best');
            end
        elseif strcmp(var, 'Inclination_deg')
            if strcmpi(fam,'Circular')
                rng = cfg.circular_inc_range;
            elseif strcmpi(fam,'Eccentric')
                rng = cfg.eccentric_inc_range;
            elseif strcmpi(fam,'Polar')
                rng = cfg.polar_inc_range;
            elseif startsWith(fam,'Frozen')
                ff = fam;
                if isfield(cfg.frozen_families, ff)
                    inc0 = cfg.frozen_families.(ff).inc_deg;
                    dt = cfg.frozen_families.(ff).inc_tolerance_deg;
                    rng = [inc0-dt, inc0+dt];
                else
                    rng = [];
                end
            else
                rng = [];
            end
            if ~isempty(rng)
                yl = ylim;
                plot([rng(1) rng(1)], yl, 'k--', 'LineWidth', 1);
                plot([rng(2) rng(2)], yl, 'k--', 'LineWidth', 1);
                ylim(yl);
                legend({'Samples','Default range'}, 'Location','best');
            end
        elseif strcmp(var, 'Eccentricity')
            if startsWith(fam,'Frozen')
                ff = fam;
                if isfield(cfg.frozen_families, ff) && isfield(cfg.frozen_families.(ff),'ecc_range')
                    rng = cfg.frozen_families.(ff).ecc_range;
                else
                    rng = [];
                end
            else
                rng = [];
            end
            if ~isempty(rng)
                yl = ylim;
                plot([rng(1) rng(1)], yl, 'k--', 'LineWidth', 1);
                plot([rng(2) rng(2)], yl, 'k--', 'LineWidth', 1);
                ylim(yl);
                legend({'Samples','Default range'}, 'Location','best');
            end
        end
    end
    sgtitle(['Population distribution: ' labels{vi}]);
    saveas(gcf, fullfile(outDir, sprintf('dist_%s_%s.png', var, ts)));
end

% Overall summary
fprintf('\nSummary by family:\n');
for fi = 1:numel(fam_list)
    fam = fam_list{fi}; sel = strcmp(all_primaries.Family, fam);
    fprintf('  %-16s N=%d\n', fam, sum(sel));
end

fprintf('\nSaved figures to %s\n', outDir);

%% Secondary plots: combined overlays and CDFs per variable
% Combined overlay histograms (all families in one axes per variable)
for vi = 1:numel(vars)
    var = vars{vi};
    figure('Name', ['Overlay Histogram: ' var], 'NumberTitle', 'off', 'Position', [100,100,1100,700]);
    hold on; grid on;
    % consistent bins across all data
    all_data = all_primaries.(var);
    if isfloat(all_data)
        q = quantile(all_data, [0.01 0.99]);
        edges = linspace(q(1), q(2), 40);
    else
        edges = []; % fallback lets histogram choose
    end
    legends = cell(1, numel(fam_list));
    for fi = 1:numel(fam_list)
        fam = fam_list{fi}; sel = strcmp(all_primaries.Family, fam);
        data = all_primaries.(var)(sel);
        if isempty(edges)
            h = histogram(data, 'Normalization','pdf', 'DisplayStyle','stairs', 'EdgeColor', colors(fi,:), 'LineWidth', 1.5);
        else
            h = histogram(data, edges, 'Normalization','pdf', 'DisplayStyle','stairs', 'EdgeColor', colors(fi,:), 'LineWidth', 1.5);
        end %#ok<NASGU>
        legends{fi} = fam;
    end
    xlabel(labels{vi}); ylabel('PDF (normalized)');
    title(['Overlay histogram by family: ' labels{vi}]);
    legend(legends, 'Location','best');
    saveas(gcf, fullfile(outDir, sprintf('overlay_hist_%s_%s.png', var, ts)));
end

% Combined overlay empirical CDFs (no toolbox dependency)
for vi = 1:numel(vars)
    var = vars{vi};
    figure('Name', ['Overlay CDF: ' var], 'NumberTitle', 'off', 'Position', [100,100,1100,700]);
    hold on; grid on;
    legends = cell(1, numel(fam_list));
    for fi = 1:numel(fam_list)
        fam = fam_list{fi}; sel = strcmp(all_primaries.Family, fam);
        data = all_primaries.(var)(sel);
        data = data(~isnan(data));
        if isempty(data)
            continue;
        end
        data = sort(data(:));
        y = (1:numel(data))' / numel(data);
        plot(data, y, 'Color', colors(fi,:), 'LineWidth', 1.5);
        legends{fi} = fam;
    end
    xlabel(labels{vi}); ylabel('Empirical CDF');
    title(['Overlay CDF by family: ' labels{vi}]);
    legend(legends, 'Location','best');
    saveas(gcf, fullfile(outDir, sprintf('overlay_cdf_%s_%s.png', var, ts)));
end

% Pairwise scatter plots across variables, colored by family
pairs = { {'Altitude_km','Inclination_deg'}, {'Altitude_km','Eccentricity'}, {'Eccentricity','Inclination_deg'} };
for pi = 1:numel(pairs)
    xvar = pairs{pi}{1}; yvar = pairs{pi}{2};
    % find labels
    xi = find(strcmp(vars, xvar)); yi = find(strcmp(vars, yvar));
    xlab = labels{xi}; ylab = labels{yi};
    figure('Name', sprintf('Scatter: %s vs %s', yvar, xvar), 'NumberTitle', 'off', 'Position', [120,120,1100,700]);
    hold on; grid on;
    legends = cell(1, numel(fam_list));
    for fi = 1:numel(fam_list)
        fam = fam_list{fi}; sel = strcmp(all_primaries.Family, fam);
        x = all_primaries.(xvar)(sel); y = all_primaries.(yvar)(sel);
        plot(x, y, '.', 'Color', colors(fi,:), 'MarkerSize', 6);
        legends{fi} = fam;
    end
    xlabel(xlab); ylabel(ylab);
    title(sprintf('Pairwise scatter by family: %s vs %s', ylab, xlab));
    legend(legends, 'Location','bestoutside');
    saveas(gcf, fullfile(outDir, sprintf('scatter_%s_vs_%s_%s.png', yvar, xvar, ts)));
end
