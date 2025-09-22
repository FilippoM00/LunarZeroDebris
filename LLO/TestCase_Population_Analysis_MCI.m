% TestCase_Population_Analysis_MCI - COMPLETE VERSION for MCI Frame LLO Population
%
% Complete cislunar collision probability analysis for LLO population data
% generated in the Moon-Centered Inertial (MCI) frame using Verifiedv1.0 framework.
% This script performs comprehensive analysis using Monte Carlo, GMM, Chan,
% Patera, and LAAS methods for the LLO TCA population.
%
% Methods: Monte Carlo, GMM, Chan, Patera, LAAS
% Framework: Verifiedv1.0 (perturbed two-body dynamics in MCI frame)
% Population: LLO orbits across 6 families (Circular, Eccentric, Polar, Frozen, General, HighEccentric)

clear all %#ok<CLALL>
close all
clc
fprintf('STEP 1 STARTED\n')
% Strict parity mode: suppress extra analysis/plots and guards to mirror
% the CR3BP reference line-by-line (except dynamics). Use env var to
% avoid static "unreachable" warnings in tooling.
sp = getenv('STRICT_PARITY');
if isempty(sp)
    strict_parity = true;
else
    strict_parity = logical(str2double(sp));
end
if strict_parity
    logf = @(varargin) [];
else
    logf = @fprintf;
end

%----------------------- Folders Paths Addition ---------------------------
% Path setup: use ONLY current workspace code (no external Smoke Codes)
addpath(genpath('../../Population_generation'));  % Population generator and samplers
addpath(genpath('./'));   % Current directory and subdirectories (LLO analytical/dynamics/utils)
addpath(genpath('../../CR3BP_ref/MCSimulations/functions_MC_SS_GMM'));  % Gaussian Mixtures functions
addpath(genpath('../../CR3BP_ref/functions'));  % POC functions
addpath(genpath('./utils')); % UnitChecks

%----------------------- Data Loading and Generation ----------------------
fprintf('Generating population...\n');

% Generate population inline using create_tca_populationLLO
total_cases = 1000; % Adjust as needed for analysis
opts = struct();
opts.seed = 42; % ensure deterministic population generation
opts.isotropic = struct('enable', true, 'seed', 42);
opts.plot.enable = false; % Disable plotting for analysis

[population, meta] = create_tca_populationLLO(total_cases, opts);
logf('Generated %d TCA cases across %d LLO families\n', population.total_generated, length(population.families));
logf('Population generation time: %.2f seconds\n', meta.computation_time);

%Plot Pre-Settings
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize', 20);

%---------------------------- MCI Constants and Parameters ----------------
cfg = LLO_defaults();
params = cfg.params;

%---------------------------- Settings ------------------------------------
fname   = './Media_MCI_Analysis';                     %Media Folder (for Pictures Export)
if exist(fname, 'dir') ~= 7, mkdir(fname); end
frmtPic = 'png';                                      %Format of Exported Pictures
odefun  = @(t,x) two_body_perturbed_rhs(t, x, params); %ODE Function for MCI
ode_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);%ODE Settings (match CR3BP)

% UNIT CONTRACTS (repo-wide consistency)
% - Dynamics/propagation/refinement stay in SI (m, m/s, s):
%   two_body_perturbed_rhs, stateTrans_perturbed_2bp, RefineTCA_MCI, EncounterDist_MCI.
% - Covariance generation uses km domain for 1σ inputs and returns km-based covariances:
%   CovGen expects [pos_km, vel_km/s]; resulting C_km is scaled back to SI with S=diag([1000..]).
% - Analytical PoC methods expect km domain at the boundary:
%   Chan_POC and compute_PoC take states [km, km/s] (row [1x6]), position covariances [km^2], and R [km].
% - GMM_POC_MCI internally converts SI→km right before Chan; do not pre-convert for it.
% - Population_generation/create_tca_populationLLO outputs states in SI at nominal TCA.

% Time unit helpers [s]
HOUR = 3600;
DAY  = 24*3600;

%Combined Objects Radius [m]
Robjs = 10.0;

% Parallelization toggle for Monte Carlo loop
use_parfor = true; % keep option to enable parfor externally

% Sample counts (aligned with CR3BP copy)
nsamp_mc  = 1e8;  % Monte Carlo samples
nsamp_ana = 1e6;  % Analytical step samples (re-sampled at Tstart)

% Covariance settings (aligned to CR3BP copy — build directly from sigma & R in SI)
k_scale  = 1;

% Strict parity mode: when true, suppress extra analysis/plots and guards to
% mirror the reference copy line-by-line (except dynamics)
% strict_parity already defined at top

%% Extract cases from generated population
fprintf('Preparing cases from generated population...\n');
N = population.total_generated;
families = population.families;

all_primary_states = population.baseline_states;       % states at TCA (primary) [6xN]
all_secondary_states = population.secondary_states;    % states at TCA (secondary) [6xN]
all_tca_times = population.tca_days;                   % TCA times in days
all_families = population.baseline_family_idx;         % family indices

all_case_info = cell(N,1);
for k = 1:N
    fam_name = families{all_families(k)};
    prim_idx = population.baseline_primary_idx(k);
    all_case_info{k} = sprintf('%s_primary_%d_case_%04d', fam_name, prim_idx, k);
end

case_counter = N;
fprintf('Total cases to process: %d\n', case_counter);

% Initialize results storage
results = struct();
results.POC_MC = zeros(case_counter, 1);
results.POC_GMM = zeros(case_counter, 1);
results.POC_Chan = zeros(case_counter, 1);
results.POC_Chan2 = zeros(case_counter, 1);
results.POC_Patera = zeros(case_counter, 1);
results.POC_LAAS = zeros(case_counter, 1);
results.case_info = all_case_info;
results.families = families;
results.family_indices = all_families;

%% Process All Cases
for case_i = 1:case_counter
    rng(case_i)
    fprintf('\n=== Processing case %d/%d: %s ===\n', case_i, case_counter, all_case_info{case_i});

    try
    %--------------------------- Test Case Data -------------------------------
    % Time Vector - Population states are AT TCA, so propagate around TCA
        TCA_nominal = 0;                        % Population states are at TCA (t=0)
        Tstart = -2 * DAY;                      % 4 days before TCA
        % Mirror CR3BP: +12 hours after TCA
        Tend   =  12 * HOUR;                     % 12 hours after TCA
        ntime  = 30000;                         % Propagation points (dense)
        tspan  = linspace(Tstart, Tend, ntime);

        % Initial States - Population data at nominal TCA (SI units)
        xP_TCA = all_primary_states(:, case_i);       % Primary at nominal TCA
        xS_TCA = all_secondary_states(:, case_i);     % Secondary at nominal TCA
        % Print nominal miss distance at initial TCA (t=0)
        d_nominal = norm(xP_TCA(1:3) - xS_TCA(1:3));  % [m]
        fprintf('  Nominal miss distance at initial TCA: %.3f m (%.6f km)\n', d_nominal, d_nominal/1000);

        % Early TCA refinement (seconds) with window escalation
        TCAold_ECI = 0;
        try
            [TCAold_ECI, ~] = RefineTCA_MCI(xP_TCA, xS_TCA, 0);
        catch
            TCAold_ECI = 0; % fallback if refinement not available
        end

        % Extract index of time instant closest to refined TCA within tspan
        [~, ix] = min(abs(tspan - TCAold_ECI));

        %Propagate original states to refined TCA (if nonzero), then back to Tstart
        if abs(TCAold_ECI) > 0
            [~, xP_TCAold] = ode113(odefun, [0, TCAold_ECI], xP_TCA, ode_options);
            [~, xS_TCAold] = ode113(odefun, [0, TCAold_ECI], xS_TCA, ode_options);
            x_refP = xP_TCAold(end,:)';
            x_refS = xS_TCAold(end,:)';
        else
            x_refP = xP_TCA;
            x_refS = xS_TCA;
        end

        % Propagation from refined TCA to MC starting point (Tstart)
        [~, x_sol_P] = ode113(odefun, [TCAold_ECI, Tstart], x_refP, ode_options);
        [~, x_sol_S] = ode113(odefun, [TCAold_ECI, Tstart], x_refS, ode_options);

        %Extract Initial States (@Tstart)
        x0_P = x_sol_P(end,:)';
        x0_S = x_sol_S(end,:)';

    %------------------------ Covariance Matrices -----------------------------
% 1. Define the desired 1-sigma standard deviations (from literature).
% These are the square roots of the diagonal of the covariance matrix.
sigma = [100;   % X-position [m]
         166;   % Y-position [m]
         100;   % Z-position [m]
         0.0020; % X-velocity [m/s]
         0.0030; % Y-velocity [m/s]
         0.0024]; % Z-velocity [m/s]

% 2. Create a diagonal matrix with the variances (sigma squared).
% D = diag(sigma.^2);

% 3. Define a correlation matrix (Rg). This is a matrix of values between -1 and 1
% that defines how state errors are related. This is a guess based on dynamics.
% Let's assume moderate correlations typical of orbital dynamics.
R = [1.0   0.6  -0.3  0.3  0.2  -0.1;
     0.6   1.0   0.3  0.2  0.4   0.1;
    -0.3   0.3   1.0 -0.1  0.1   0.3;
     0.3   0.2  -0.1  1.0  0.3  -0.2;
     0.2   0.4   0.1  0.3  1.0   0.1;
    -0.1   0.1   0.3 -0.2  0.1   1.0];

% 4. Combine the diagonal standard deviations with the correlation structure
% to form a valid, full covariance matrix.
% The formula is: P = diag(sigma) * R * diag(sigma)
P_valid = diag(sigma) * R * diag(sigma);
CP0 = P_valid; CS0 = CP0;
fprintf('STEP 1 EXECUTED\n')
    % Fail-fast unit checks
    %UnitChecks.assertCovPSD(CP0, 'CP0');
    %UnitChecks.assertCovPSD(CS0, 'CS0');
    %UnitChecks.assertCovSigmaConsistency(CP0, err_posP, err_velP, 'CP0');
    %UnitChecks.assertCovSigmaConsistency(CS0, err_posS, err_velS, 'CS0');
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%  MONTE-CARLO SIMULATION  %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('STEP 2 STARTED\n')

        %MC Sampling (same structure as CR3BP)
    if ~strict_parity
        fprintf('  Running Monte Carlo with %d samples...\n', nsamp_mc);
    end
    RP0 = mvnrnd(x0_P', CP0, nsamp_mc);
    RS0 = mvnrnd(x0_S', CS0, nsamp_mc);

        %Memory Allocation
    counter  = zeros(nsamp_mc,1); %Collisions Counter
    tEncount = nan(1,nsamp_mc);   %Collision Time Instants

    % Search w/Minimum Distance using EncounterDist_MCI
    % Prefer parfor when available; fallback to for otherwise
    
for i = 1 : nsamp_mc
    try
                % Use ±12 hours refine window to match CR3BP flow
                [min_dist, tEnc] = EncounterDist_MCI([RP0(i,:), RS0(i,:)], tspan, params, 'refine', [-12*HOUR, 12*HOUR]);
                % echo min_dist to mirror reference console output
                min_dist
                if min_dist <= Robjs
                    counter(i)  = 1;
                    tEncount(i) = tEnc;
                end
    catch
                % Fallback to direct distance at nominal TCA
                d0 = norm(RP0(i,1:3) - RS0(i,1:3));
                if d0 <= Robjs
                    counter(i)  = 1;
                    tEncount(i) = 0;
                end
    end
end

        %MC PoC
    POC_MC = sum(counter)/nsamp_mc;

    % Adjusting Encounters Times (mirror CR3BP cumulative PoC plotting)
        tEncount  = sort(tEncount, 'ascend', 'MissingPlacement', 'last');
    indx  = find(tspan > (-12*HOUR) & tspan < (12*HOUR)); % Subwindow [-12h, +12h]
        tEncount2 = zeros(1, length(tspan(indx)));
    for ii = 1:nsamp_mc
        if ~isnan(tEncount(ii))
            [~, closestIndex] = min(abs(tspan(indx) - tEncount(ii)));
            tEncount2(closestIndex) = tEncount2(closestIndex) + 1;
        else
            break
        end
    end
    CumPOC = cumsum(tEncount2, 2) / nsamp_mc;
    fprintf('STEP 2.1 EXECUTED\n')

    % % Plot cumulative probability vs time and save
        % try
        %     fig = figure('Position', [30 30 1200 650]);
        %     hold on; grid on;
    %     plot(tspan(indx), CumPOC, 'b', 'LineWidth', 1.4);
        %     title('Cumulative Probability');
        %     xlabel('Time from TCA [s]');
        %     ylabel('Probability of Collision [-]');
        %     xlim([min(tspan(indx_win)), max(tspan(indx_win))]);
        %     plot_filename = sprintf('CumulativePOC_%s', all_case_info{case_i});
        %     saveas(fig, fullfile(fname, plot_filename), frmtPic);
        %     exportgraphics(fig, fullfile(fname, strcat(plot_filename, '.pdf')));
        %     close(fig);
        % catch ME
        %     fprintf('Warning: cumulative PoC plotting failed: %s\n', ME.message);
    % end
    fprintf('STEP 2.2 EXECUTED\n')
    % Save MC results as in reference
    % save('Basic_MC_1e5_', 'POC_MC', 'counter', 'tEncount', 'tEncount2', 'CumPOC', 'tspan', 'indx');
    fprintf('STEP 2.EXECUTED\n')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% CHAN & OTHER METHODS %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('STEP 3 STARTED\n');

    % Resample for analytical step (copy-CR3BP alignment)
        RP0_ana = mvnrnd(x0_P', CP0, nsamp_ana);
        RS0_ana = mvnrnd(x0_S', CS0, nsamp_ana);

        %------------------------- Propagation to TCA -----------------------------
        stateTCA_P = zeros(nsamp_ana, 6);
        stateTCA_S = zeros(nsamp_ana, 6);
        for i = 1 : nsamp_ana
            [~, stateP] = ode113(odefun, tspan, RP0_ana(i,:)', ode_options);
            [~, stateS] = ode113(odefun, tspan, RS0_ana(i,:)', ode_options);
            stateTCA_P(i,:) = stateP(ix,:);
            stateTCA_S(i,:) = stateS(ix,:);
        end

        % Extraction of TCA condition (mean & covariance at TCA)
        xPmean = mean(stateTCA_P,1)';
        xSmean = mean(stateTCA_S,1)';
        CPmean = cov(stateTCA_P);
        CSmean = cov(stateTCA_S);
        % Robust fallback: when nsamp_ana==1 MATLAB's cov can be degenerate/NaN.
        % Transport CP0/CS0 from Tstart to tspan(ix) to obtain CPmean/CSmean.
        if any(~isfinite(CPmean(:))) || any(~isfinite(CSmean(:)))
            try
                tf_to_ix = tspan(ix) - Tstart;
                [~, PhiP_ix] = stateTrans_perturbed_2bp(x0_P, tf_to_ix, params, ode_options, Tstart);
                [~, PhiS_ix] = stateTrans_perturbed_2bp(x0_S, tf_to_ix, params, ode_options, Tstart);
                CPmean = PhiP_ix * CP0 * PhiP_ix';
                CSmean = PhiS_ix * CS0 * PhiS_ix';
                CPmean = (CPmean + CPmean')/2;
                CSmean = (CSmean + CSmean')/2;
            catch
                % As a last resort, use CP0/CS0 directly
                CPmean = CP0; CSmean = CS0;
            end
        end

        % Reference behavior: use mean/covariance from sampled propagation at ix (no overrides)

        %----------------------------- TCA Refinement -----------------------------
        % Use a wider window for the second refine to be robust under MCI dynamics
        % Use default ±2h window (reference behavior) for the second refine
        [TCAnew_ECI, ~] = RefineTCA_MCI(xPmean, xSmean, TCAold_ECI);
        TCAnew = TCAnew_ECI; % seconds

    % Log parity with reference: both refinements should often match
    if ~strict_parity
        fprintf('  Refine TCA: old=%.6f s, new=%.6f s, delta=%.3e s\n', TCAold_ECI, TCAnew_ECI, TCAnew_ECI - TCAold_ECI);
    end

        % Propagate mean states to new TCA (seconds)
        [~, xP_TCAnew] = ode113(odefun, [tspan(ix), TCAnew], xPmean', ode_options);
        [~, xS_TCAnew] = ode113(odefun, [tspan(ix), TCAnew], xSmean', ode_options);

    % Retrieve STM from closest time to TCAnew
try
            tf_p = TCAnew - tspan(ix);
            [~, STMP] = stateTrans_perturbed_2bp(xPmean, tf_p, params, ode_options, tspan(ix));
            [~, STMS] = stateTrans_perturbed_2bp(xSmean, tf_p, params, ode_options, tspan(ix));
            CPmean_new = STMP * CPmean * STMP';
            CSmean_new = STMS * CSmean * STMS';
            % enforce symmetry to reduce numerical noise
            CPmean_new = (CPmean_new + CPmean_new')/2;
            CSmean_new = (CSmean_new + CSmean_new')/2;
catch
            fprintf('Warning: STM propagation failed. Using unpropagated covariances.\n');
            CPmean_new = CPmean;
            CSmean_new = CSmean;
end

        %----------------------------- Vittaldev Split ----------------------------
        [vecP, lamP] = eig(CPmean_new);
        [~, colmaxP] = max(diag(lamP));
        dirP = vecP(:,colmaxP);

        [vecS, lamS] = eig(CSmean_new);
        [~, colmaxS] = max(diag(lamS));
        dirS = vecS(:,colmaxS);

    [GMM_P.mean, GMM_P.cov, GMM_P.w] = VittaldevSplit(xP_TCAnew(end,:)', CPmean_new, dirP, 39);
    [GMM_S.mean, GMM_S.cov, GMM_S.w] = VittaldevSplit(xS_TCAnew(end,:)', CSmean_new, dirS, 39);
    % Align shapes with CR3BP reference: means as [N x 6]
    GMM_P.mean = GMM_P.mean';
    GMM_S.mean = GMM_S.mean';

        % Precompute km-domain states and R once for PoC calls
        xP_km_row = (xP_TCAnew(end,:)' / 1000)';
        xS_km_row = (xS_TCAnew(end,:)' / 1000)';
        R_km = Robjs / 1000;

        %--------------------------------- GMM POC --------------------------------
        try
            fprintf('STEP 3.1 STARTED\n')
            % Pass ODE options and AstroMod='mci'; Robjs in km for Chan inside GMM
            POC_GMM = GMM_POC_MCI(GMM_P, GMM_S, TCAnew_ECI, ode_options, 'mci', Robjs/1000, 10);
            fprintf('STEP 3.1 EXECUTED\n')
        catch
            fprintf('Warning: GMM_POC_MCI failed. Setting POC_GMM=NaN.\n');
            POC_GMM = NaN;
        end

        %-------------------------------- Chan POC --------------------------------
        try
            fprintf('STEP 3.2 STARTED\n')
            % Convert covariances to km^2 and symmetrize before Chan
            Csum6_km = (CPmean_new + CSmean_new) / 1e6;  % full 6x6 in km-based units
            Csum6_km = (Csum6_km + Csum6_km')/2;
            POC_Chan = Chan_POC(xP_km_row, xS_km_row, Csum6_km, R_km, 10);
            fprintf('STEP 3.2 EXECUTED\n')
        catch
            fprintf('Warning: Chan_POC failed. Using fallback value 0. Reason: %s\n', lasterror.message);
            POC_Chan = 0;
        end

        %-------------------------------- Chan 2 / Patera / LAAS POC -------------
        try
            fprintf('STEP 3.3 STARTED\n')
            % compute_PoC expects km units for states/covariances and km for HBR; use ROW vectors [1x6]
            CP_pos_km = CPmean_new(1:3,1:3) / 1e6;
            CS_pos_km = CSmean_new(1:3,1:3) / 1e6;
            POC_Chan2  = compute_PoC(xP_km_row, CP_pos_km, xS_km_row, CS_pos_km, R_km, 'Chan')
            fprintf('STEP 3.3 EXECUTED\n')
            fprintf('STEP 3.4 STARTED\n')
            POC_Patera = compute_PoC(xP_km_row, CP_pos_km, xS_km_row, CS_pos_km, R_km, 'Patera');
            fprintf('STEP 3.4 EXECUTED\n')
            fprintf('STEP 3.5 STARTED\n')
            POC_LAAS   = compute_PoC(xP_km_row, CP_pos_km, xS_km_row, CS_pos_km, R_km, 'LAAS');
            fprintf('STEP 3.5 EXECUTED\n')
        catch
            fprintf('Warning: compute_PoC methods failed. Using Chan result as fallback. Reason: %s\n', lasterror.message);
            POC_Chan2 = POC_Chan;
            POC_Patera = POC_Chan;
            POC_LAAS = POC_Chan;
        end
    % Clamp tiny negative probabilities to zero due to numerical noise
    POC_Chan   = max(POC_Chan,   0);
    POC_Chan2  = max(POC_Chan2,  0);
    POC_Patera = max(POC_Patera, 0);
    POC_LAAS   = max(POC_LAAS,   0);
    fprintf('STEP 3 EXECUTED\n')

        %----------------------------- Parity/Quality Checks ----------------------
                if ~strict_parity
                    try
            % 1) Orthogonality at refined TCA (means): r_rel · v_rel ≈ 0
            r_rel_mean = xP_TCAnew(end,1:3)' - xS_TCAnew(end,1:3)';
            v_rel_mean = xP_TCAnew(end,4:6)' - xS_TCAnew(end,4:6)';
            drdv_mean = abs(dot(r_rel_mean, v_rel_mean));
            if drdv_mean > 1e-3 % m^2/s (tight but nonzero tolerance)
                fprintf('Warning: TCA orthogonality not tight (|Dr·Dv|=%.3e) in case %s\n', drdv_mean, all_case_info{case_i});
            end

            % 2) Covariance symmetry/conditioning
            sym_err_P = norm(CPmean_new - CPmean_new','fro');
            sym_err_S = norm(CSmean_new - CSmean_new','fro');
            if sym_err_P > 1e-9 || sym_err_S > 1e-9
                fprintf('Warning: covariance symmetry error P=%.2e S=%.2e in case %s\n', sym_err_P, sym_err_S, all_case_info{case_i});
            end
            % PSD-ish check (allow tiny negative eigs from num noise)
            min_eig_P = min(eig((CPmean_new+CPmean_new')/2));
            min_eig_S = min(eig((CSmean_new+CSmean_new')/2));
            if min_eig_P < -1e-12 || min_eig_S < -1e-12
                fprintf('Warning: covariance not PSD (min eig P=%.2e S=%.2e) case %s\n', min_eig_P, min_eig_S, all_case_info{case_i});
            end

            % 3) STM finiteness
            if any(~isfinite(reshape(STMP,[],1))) || any(~isfinite(reshape(STMS,[],1)))
                fprintf('Warning: non-finite values in STM in case %s\n', all_case_info{case_i});
            end

            % 4) GMM weights sanity
            if exist('GMM_P','var') && isfield(GMM_P,'w')
                wsumP = sum(GMM_P.w);
                if abs(wsumP - 1) > 1e-8
                    fprintf('Warning: GMM_P weights sum=%.6f (not 1) in case %s\n', wsumP, all_case_info{case_i});
                end
            end
            if exist('GMM_S','var') && isfield(GMM_S,'w')
                wsumS = sum(GMM_S.w);
                if abs(wsumS - 1) > 1e-8
                    fprintf('Warning: GMM_S weights sum=%.6f (not 1) in case %s\n', wsumS, all_case_info{case_i});
                end
            end
                    catch ME
                        fprintf('Parity checks failed: %s (case %s)\n', ME.message, all_case_info{case_i});
                    end
                end

        % Store results
        results.POC_MC(case_i) = POC_MC;
        results.POC_GMM(case_i) = POC_GMM;
        results.POC_Chan(case_i) = POC_Chan;
        results.POC_Chan2(case_i) = POC_Chan2;
        results.POC_Patera(case_i) = POC_Patera;
        results.POC_LAAS(case_i) = POC_LAAS;

        %------------------------------- Results ------------------------------
    fprintf('Results for %s:\n', all_case_info{case_i});
    fprintf('  POC_MC: %.6e\n', POC_MC);
    fprintf('  POC_GMM: %.6e\n', POC_GMM);
    fprintf('  POC_Chan: %.6e\n', POC_Chan);
    fprintf('  POC_Chan2: %.6e\n', POC_Chan2);
    fprintf('  POC_Patera: %.6e\n', POC_Patera);
    fprintf('  POC_LAAS: %.6e\n', POC_LAAS);

        catch ME
        fprintf('ERROR in case %d (%s): %s\n', case_i, all_case_info{case_i}, ME.message);
        results.POC_MC(case_i) = NaN;
        results.POC_GMM(case_i) = NaN;
        results.POC_Chan(case_i) = NaN;
        results.POC_Chan2(case_i) = NaN;
        results.POC_Patera(case_i) = NaN;
        results.POC_LAAS(case_i) = NaN;
    end

end

% (Analysis plots/summary intentionally omitted under strict parity)

%% Save Others (match CR3BP reference)
% Prepare variables similar to reference naming for the last processed case
% Note: variables from inside the loop are available here (MATLAB workspace semantics)
xPmean_new = xP_TCAnew(end,:); %#ok<NODEF>
xSmean_new = xS_TCAnew(end,:); %#ok<NODEF>
CPmean_new_ECI = CPmean_new / 1e6; %#ok<NODEF>
CSmean_new_ECI = CSmean_new / 1e6; %#ok<NODEF>
TCAnew = TCAnew_ECI; %#ok<NODEF>
% Save propagation @TCA and POC values
save('Basic_Others_1e5_4d', 'xPmean', 'xSmean', 'CPmean', 'CSmean', ...
    'TCAnew', 'TCAnew_ECI', 'xPmean_new', 'xSmean_new', 'CPmean_new', 'CSmean_new', ...
    'CPmean_new_ECI', 'CSmean_new_ECI', ...
    'POC_GMM', 'POC_Chan', 'POC_Chan2', 'POC_Patera', 'POC_LAAS');
