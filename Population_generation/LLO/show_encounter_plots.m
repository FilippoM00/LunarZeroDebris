% SHOW_ENCOUNTER_PLOTS
% Generate a TCA population and SHOW the key encounter plots:
% - Miss distance
% - Approach angles
% - Relative velocities
% Plots are interactive and also saved under Media_MCI_Analysis.

clear; clc; close all;

% Ensure workspace paths
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
popRoot  = fileparts(thisDir); % ..\Population_generation
addpath(genpath(popRoot));

% Make sure figures are visible
set(0, 'DefaultFigureVisible', 'on');

% Output dir
outDir = fullfile(thisDir, 'Media_MCI_Analysis');
if ~exist(outDir, 'dir'); mkdir(outDir); end
ts = datestr(now,'yyyy-mm-dd_HHMMSS');

% Configure generation (disable internal plotting to avoid extra figures)
opts = struct();
opts.plot = struct('enable', false, ...
                   'saveDir', outDir, ...
                   'format', 'png');

% Total encounters to generate (adjust as needed)
total = 500;

fprintf('Generating %d encounters and showing encounter plots...\n', total);
[population, meta] = create_tca_populationLLO(total, opts); %#ok<ASGLU>

% Convenience variables
dmiss_km = population.dmiss_km(:);
vrel_kms = population.vrel_kms(:);
angles_deg = population.angles_deg(:);

%% 1) Miss distance plots
figure('Name','Miss Distance Analysis','NumberTitle','off','Position',[100,100,900,650]);
% Log histogram
subplot(2,2,1);
histogram(log10(dmiss_km), 40, 'Normalization','pdf'); grid on;
xlabel('log_{10}(Miss Distance) [km]'); ylabel('PDF');
title('Miss Distance (Log Scale)');
% Linear histogram
subplot(2,2,2);
histogram(dmiss_km, 40, 'Normalization','pdf'); grid on;
xlabel('Miss Distance [km]'); ylabel('PDF');
title('Miss Distance (Linear Scale)');
% CDF
subplot(2,2,3);
[f,x] = ecdf(dmiss_km); plot(x,f,'LineWidth',2); grid on;
xlabel('Miss Distance [km]'); ylabel('CDF'); title('Miss Distance CDF');
% Stats
subplot(2,2,4); axis off;
txt = sprintf(['Miss Distance Stats\n' ...
               'Min: %.2e km\nMax: %.2e km\nMean: %.2e km\nMedian: %.2e km\nStd: %.2e km\nN: %d'], ...
               min(dmiss_km), max(dmiss_km), mean(dmiss_km), median(dmiss_km), std(dmiss_km), numel(dmiss_km));
text(0.05,0.95,txt,'VerticalAlignment','top');
saveas(gcf, fullfile(outDir, sprintf('miss_distance_analysis_SHOWN_%s.png', ts)));

%% 2) Relative velocity plots
figure('Name','Relative Velocity Analysis','NumberTitle','off','Position',[120,120,900,650]);
% Linear histogram
subplot(2,2,1);
histogram(vrel_kms, 40, 'Normalization','pdf'); grid on;
xlabel('Relative Velocity [km/s]'); ylabel('PDF');
title('Relative Velocity');
% Log histogram (avoid -Inf)
subplot(2,2,2);
log_v = log10(max(vrel_kms, 1e-10));
histogram(log_v, 40, 'Normalization','pdf'); grid on;
xlabel('log_{10}(Relative Velocity) [km/s]'); ylabel('PDF');
title('Relative Velocity (Log Scale)');
% CDF
subplot(2,2,3);
[f,x] = ecdf(vrel_kms); plot(x,f,'LineWidth',2); grid on;
xlabel('Relative Velocity [km/s]'); ylabel('CDF'); title('Relative Velocity CDF');
% Stats
subplot(2,2,4); axis off;
txt = sprintf(['Relative Velocity Stats\n' ...
               'Min: %.4f km/s\nMax: %.4f km/s\nMean: %.4f km/s\nMedian: %.4f km/s\nStd: %.4f km/s\nN: %d'], ...
               min(vrel_kms), max(vrel_kms), mean(vrel_kms), median(vrel_kms), std(vrel_kms), numel(vrel_kms));
text(0.05,0.95,txt,'VerticalAlignment','top');
saveas(gcf, fullfile(outDir, sprintf('relative_velocity_analysis_SHOWN_%s.png', ts)));

%% 3) Encounter angles plots
figure('Name','Encounter Angles Analysis','NumberTitle','off','Position',[140,140,900,650]);
% Angular histogram
subplot(2,2,1);
histogram(angles_deg, 36, 'Normalization','pdf'); grid on;
xlabel('Encounter Angle [deg]'); ylabel('PDF');
title('Encounter Angles');
% Polar histogram
subplot(2,2,2);
polarhistogram(deg2rad(angles_deg), 36, 'Normalization','pdf');
title('Encounter Angles (Polar)');
% CDF
subplot(2,2,3);
[f,x] = ecdf(angles_deg); plot(x,f,'LineWidth',2); grid on;
xlabel('Encounter Angle [deg]'); ylabel('CDF'); title('Encounter Angles CDF');
% Stats
subplot(2,2,4); axis off;
txt = sprintf(['Encounter Angles Stats\n' ...
               'Min: %.1f°\nMax: %.1f°\nMean: %.1f°\nMedian: %.1f°\nStd: %.1f°\nN: %d'], ...
               min(angles_deg), max(angles_deg), mean(angles_deg), median(angles_deg), std(angles_deg), numel(angles_deg));
text(0.05,0.95,txt,'VerticalAlignment','top');
saveas(gcf, fullfile(outDir, sprintf('encounter_angles_analysis_SHOWN_%s.png', ts)));

fprintf('Done. Plots shown and saved to %s\n', outDir);
