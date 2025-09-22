% Test_PoC_MCI.m - Mini PoC Analysis for MCI Frame LLO Population
%
% Mirrors the full TestCase_Population_Analysis_MCI.m but with fewer cases
% and smaller MC samples for verification. Tests primary-secondary pairs
% and PoC computation.

clear all
close all
clc

%----------------------- Folders Paths Addition ---------------------------%
addpath(genpath('../../Population_generation'));
addpath(genpath('./'));
addpath(genpath('../../CR3BP_ref/MCSimulations/functions_MC_SS_GMM'));

%----------------------- Data Loading and Generation ----------------------%
fprintf('Generating small LLO TCA population for testing...\n');

% Generate small population (2-3 cases)
total_cases = 6;
opts = struct();
opts.isotropic = struct('enable', true, 'seed', 42);
opts.plot.enable = false;

[population, meta] = create_tca_populationLLO(total_cases, opts);

fprintf('Generated %d TCA cases\n', population.total_generated);

%---------------------------- MCI Constants and Parameters ----------------%
cfg = LLO_defaults();
params = cfg.params;

%---------------------------- Settings ------------------------------------%
fname = './Media_Test'; % Test media folder
if exist(fname, 'dir') ~= 7, mkdir(fname); end
frmtPic = 'png';
odefun = @(t,x) two_body_perturbed_rhs(t, x, params);
ode_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% Time unit helpers
HOUR = 3600;
DAY = 24*3600;

% Combined Objects Radius
Robjs = 10.0;

% Covariance settings (define in SI, convert for CovGen which expects km / km/s)
k_scale = 1;
err_posP = (1000.0/3) * k_scale; % [m] 1σ
err_velP = (0.01/3) * k_scale;    % [m/s] 1σ
err_posS = (1000.0/3) * k_scale;  % [m] 1σ
err_velS = (0.1/3) * k_scale;     % [m/s] 1σ

% Prepare km inputs for CovGen
err_posP_km   = err_posP / 1000;    % [km]
err_velP_kmps = err_velP / 1000;    % [km/s]
err_posS_km   = err_posS / 1000;    % [km]
err_velS_kmps = err_velS / 1000;    % [km/s]

%% Extract cases
fprintf('Preparing test cases...\n');
N = population.total_generated;
all_primary_states = population.baseline_states;
all_secondary_states = population.secondary_states;
all_tca_times = population.tca_days;
all_families = population.baseline_family_idx;

all_case_info = cell(N,1);
for k = 1:N
    fam_name = population.families{all_families(k)};
    prim_idx = population.baseline_primary_idx(k);
    all_case_info{k} = sprintf('%s_primary_%d_case_%04d', fam_name, prim_idx, k);
end

case_counter = N;
fprintf('Total test cases: %d\n', case_counter);

% Initialize results
results.POC_MC = zeros(case_counter, 1);
results.POC_Chan = zeros(case_counter, 1);

%% Process Test Cases
for case_i = 1:case_counter

    fprintf('\n=== Processing test case %d/%d: %s ===\n', case_i, case_counter, all_case_info{case_i});

    % Time Vector
    TCA_nominal = 0;
    Tstart = -4 * DAY;
    Tend = 2 * HOUR;
    ntime = 1000; % Fewer points for test
    tspan = linspace(Tstart, Tend, ntime);

    % Initial States
    xP_TCA = all_primary_states(:, case_i);
    xS_TCA = all_secondary_states(:, case_i);

    % % TCA Refinement
    % try
    %     win = 2*HOUR;
    %     [TCAold_ECI, ~] = RefineTCA_MCI(xP_TCA, xS_TCA, 0, win, win);
    % catch
         TCAold_ECI = 0;
    % end

    % Extract TCA index
    [~, ix] = min(abs(tspan - TCAold_ECI));

    % Propagate to Tstart
    % if abs(TCAold_ECI) > 0
    %     [~, xP_TCAold] = ode113(odefun, [0, TCAold_ECI], xP_TCA, ode_options);
    %     [~, xS_TCAold] = ode113(odefun, [0, TCAold_ECI], xS_TCA, ode_options);
    %     x_refP = xP_TCAold(end,:)';
    %     x_refS = xS_TCAold(end,:)';
    % else
        x_refP = xP_TCA;
        x_refS = xS_TCA;
    % end

    [~, x_sol_P] = ode113(odefun, [TCAold_ECI, Tstart], x_refP, ode_options);
    [~, x_sol_S] = ode113(odefun, [TCAold_ECI, Tstart], x_refS, ode_options);

    x0_P = x_sol_P(end,:)';
    x0_S = x_sol_S(end,:)';

    % Covariances: use CovGen (km units) and convert back to SI (m, m/s)
    [CP0_km, ~] = CovGen([err_posP_km, err_velP_kmps], 30041993 + case_i);
    [CS0_km, ~] = CovGen([err_posS_km, err_velS_kmps], 28081999 + case_i);
    S = diag([1000,1000,1000, 1000,1000,1000]);
    CP0 = S * CP0_km * S';
    CS0 = S * CS0_km * S';

    %%%%%%%%%%%%%%%%%%%%%%%%%  MONTE-CARLO SIMULATION  %%%%%%%%%%%%%%%%%%%%%%%%
    nsamp = 1000; % Small sample for test
    fprintf('  Running Monte Carlo with %d samples...\n', nsamp);
    RP0 = mvnrnd(x0_P', CP0, nsamp);
    RS0 = mvnrnd(x0_S', CS0, nsamp);

    counter = zeros(nsamp,1);
    tEncount = nan(1,nsamp);

    for i = 1:nsamp
        try
            [min_dist, tEnc] = EncounterDist_MCI([RP0(i,:), RS0(i,:)], tspan, params, 'refine', [-2*HOUR, 2*HOUR]);
            min_dist
            tEnc
            if min_dist <= Robjs
                counter(i) = 1;
                tEncount(i) = tEnc;
            end
        catch
            d0 = norm(RP0(i,1:3) - RS0(i,1:3));
            if d0 <= Robjs
                counter(i) = 1;
                tEncount(i) = 0;
            end
        end
    end

    POC_MC = sum(counter)/nsamp;

    % Cumulative PoC plot
    tEncount = sort(tEncount, 'ascend', 'MissingPlacement', 'last');
    indx_win = find(tspan > (-2*HOUR) & tspan < (2*HOUR));
    tEncount2 = zeros(1, length(tspan(indx_win)));
    for ii = 1:nsamp
        if ~isnan(tEncount(ii))
            [~, closestIndex] = min(abs(tspan(indx_win) - tEncount(ii)));
            tEncount2(closestIndex) = tEncount2(closestIndex) + 1;
        else
            break
        end
    end
    CumPOC = cumsum(tEncount2, 2) / nsamp;

    fig = figure('Position', [30 30 800 400]);
    plot(tspan(indx_win)/3600, CumPOC, 'b', 'LineWidth', 1.4);
    title(sprintf('Cumulative PoC - %s', all_case_info{case_i}));
    xlabel('Time from TCA [hours]');
    ylabel('Probability of Collision [-]');
    xlim([min(tspan(indx_win))/3600, max(tspan(indx_win))/3600]);
    grid on;
    saveas(fig, fullfile(fname, sprintf('CumPOC_Test_%s.png', all_case_info{case_i})));
    close(fig);

    %%%%%%%%%%%%%%%%%%%%%%%%% CHAN METHOD %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('  Running Chan method...\n');

    % Propagation to TCA
    stateTCA_P = zeros(nsamp, 6);
    stateTCA_S = zeros(nsamp, 6);

    for i = 1:nsamp
        [~, stateP] = ode113(odefun, tspan, RP0(i,:)', ode_options);
        [~, stateS] = ode113(odefun, tspan, RS0(i,:)', ode_options);
        stateTCA_P(i,:) = stateP(ix,:);
        stateTCA_S(i,:) = stateS(ix,:);
    end

    xPmean = mean(stateTCA_P,1)';
    xSmean = mean(stateTCA_S,1)';
    CPmean = cov(stateTCA_P);
    CSmean = cov(stateTCA_S);

    % TCA Refinement
    [TCAnew_ECI, ~] = RefineTCA_MCI(xPmean, xSmean, TCAold_ECI, 2*HOUR, 2*HOUR);
    TCAnew = TCAnew_ECI;

    % Propagate means
    [~, xP_TCAnew] = ode113(odefun, [tspan(ix), TCAnew], xPmean', ode_options);
    [~, xS_TCAnew] = ode113(odefun, [tspan(ix), TCAnew], xSmean', ode_options);

    % STM for covariances
    try
        tf_p = TCAnew - tspan(ix);
        [~, STMP] = stateTrans_perturbed_2bp(xPmean, tf_p, params, ode_options, tspan(ix));
        [~, STMS] = stateTrans_perturbed_2bp(xSmean, tf_p, params, ode_options, tspan(ix));
        CPmean_new = STMP * CPmean * STMP';
        CSmean_new = STMS * CSmean * STMS';
    catch
        CPmean_new = CPmean;
        CSmean_new = CSmean;
    end

    % Chan PoC
    try
        xP_km = (xP_TCAnew(end,:)' / 1000)';
        xS_km = (xS_TCAnew(end,:)' / 1000)';
        C_km = (CPmean_new(1:3,1:3) + CSmean_new(1:3,1:3)) / 1e6;
        R_km = Robjs / 1000;
        POC_Chan = Chan_POC(xP_km, xS_km, C_km, R_km, 10);
    catch
        POC_Chan = 0;
    end

    % Store results
    results.POC_MC(case_i) = POC_MC;
    results.POC_Chan(case_i) = POC_Chan;

    fprintf('  Results: POC_MC = %.6e, POC_Chan = %.6e\n', POC_MC, POC_Chan);

end

%% Summary
fprintf('\n=== TEST SUMMARY ===\n');
fprintf('Tested %d cases with %d MC samples each\n', case_counter, nsamp);
fprintf('PoC Results:\n');
for i = 1:case_counter
    fprintf('  %s: MC=%.6e, Chan=%.6e\n', all_case_info{i}, results.POC_MC(i), results.POC_Chan(i));
end

fprintf('Test complete. Check plots in %s\n', fname);
