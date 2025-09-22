% Smoke test: create a small TCA population for LLO primaries and plot basics

clearvars; close all; clc
addpath(genpath('.'));

% Small run
N = 200;
opts = struct();
opts.seed = 12345;
opts.time = struct('min_days', 1.0, 'max_days', 4.0);
opts.miss = struct('min_km', 1e-6, 'max_km', 3.0);
opts.isotropic = struct('enable', true, 'seed', 7);
opts.units = struct('LU_km', 1.0, 'TU_days', 1.0/86400.0); % 1 km, 1 sec

[pop, meta] = create_tca_population_LLO(N, opts); %#ok<NASGU>

% Simple sanity plots
figure('Name','LLO TCA: vrel & angles','Position',[100,100,1000,420]);
subplot(1,2,1); histogram(pop.vrel_kms, 40); grid on; xlabel('v_{rel} [km/s]'); ylabel('count');
subplot(1,2,2); histogram(pop.angles_deg, 0:5:360); grid on; xlim([0,360]); xlabel('approach angle [deg]'); ylabel('count');

% Check orthogonality numerically
dd = arrayfun(@(k) dot(pop.r_rel(:,k), pop.v_rel(:,k)), 1:pop.n_encounters);
fprintf('Max |r_rel·v_rel| = %.3e (LU^2/TU)\n', max(abs(dd)));

% Save a quick figure
outdir = fullfile('.', 'Media'); if ~exist(outdir,'dir'), mkdir(outdir); end
saveas(gcf, fullfile(outdir, 'smoke_LLO_TCA_vrel_angles.png'));
fprintf('Saved plot to %s\n', fullfile(outdir, 'smoke_LLO_TCA_vrel_angles.png'));
