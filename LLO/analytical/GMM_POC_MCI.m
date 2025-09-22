function [POC_GMM] = GMM_POC_MCI(GMM_p, GMM_s, t0, options, AstroMod, Robjs, ChanOrd)
%
% Probability of Collision w/Gaussian Mixtures (MCI variant)
%
%DESCRIPTION:
%This code provides the implementation of Gaussian Mixtures Method for the
%computation of POC for both Keplerian and MCI (inertial) cases.
%
%PROTOTYPE
%   [POC_GMM] = GMM_POC_MCI(GMM_p, GMM_s, t0, options, AstroMod, Robjs, ChanOrd)
%
%--------------------------------------------------------------------------
% INPUTS:
%   GMM_p      [1x1]       Primary GMM                [-] (struct)
%   GMM_s      [1x1]       Secondary GMM              [-] (struct)
%   t0         [1x1]       Nominal TCA                [s]
%   options    [---]       ODE Options                [-]
%   AstroMod   [1x1]       Astrodynamical Model       [string] ('kep' or 'mci')
%   Robjs      [1x1]       Combined Objects Radii     [km] (or [-])
%   ChanOrd    [1x1]       Chan Order of Expansion    [-]
%--------------------------------------------------------------------------
% OUTPUTS:
%   POC_GMM    [1x1]       Probability of Collision w/GMM  [-]
%--------------------------------------------------------------------------
%
%NOTES:
% - The input "t0" shall be the nominal TCA of the distribution. It is used
%   to refine the TCA of the mixtures couples (see DrDv_MCI.m / RefineTCA_MCI.m).
%
%CALLED FUNCTIONS:
% DrDv, DrDv_MCI, analytic_kep_prop, Chan_POC, RefineTCA, RefineTCA_MCI,
% stateTransKeplerian, stateTrans_perturbed_2bp
%
%AUTHOR(s):
%Luigi De Maria, 2022 (adapted to MCI)
%

%% Main Code

%POC Initialization
POC_GMM = 0;

% Small optimization: prepare MCI params and ODE once (no-op for 'kep')
if strcmpi(AstroMod,'mci')
    cfg = LLO_defaults();
    params = cfg.params;
    ode_fun = @(t,x) two_body_perturbed_rhs(t, x, params);
end

% Pre-fetch weights (no pruning; keep exact parity with reference)
wP = GMM_p.w(:);
wS = GMM_s.w(:);

%Looping Over the GMEs Couples to Compute the Overall POC w/GMM
for i = 1:length(wP)
    for j = 1:length(wS)
        %Refine Mixtures Couple TCA
        switch AstroMod
            %------------------------- KEPLERIAN --------------------------
            case 'kep'
                %Constants
                mu   = 4902.800076;  %Moon Planetary Parameter [km3][s-2] (MODIFIED for lunar analysis)
                %Definition of Cost Function
                fun  = @(t) DrDv(t0, t, GMM_p.mean(i,:), GMM_s.mean(j,:), options, 'kep');
                %fsolve Options
                fsolveopt = optimoptions('fsolve','FunctionTolerance',1e-30,'Algorithm',...
                                         'levenberg-marquardt','Display','off');
                %Real TCA
                tTCA = fsolve(fun, t0+0.1, fsolveopt);
                %Propagation to Refined TCA of States
                [xP] = analytic_kep_prop(GMM_p.mean(i,:)', [t0 tTCA], mu); xP = xP';
                [xS] = analytic_kep_prop(GMM_s.mean(j,:)', [t0 tTCA], mu); xS = xS';
                
                %Propagation to Refined TCA of Covariance Matrix
                STMP = stateTransKeplerian([t0 tTCA], GMM_p.mean(i,:)', mu, options);
                STMS = stateTransKeplerian([t0 tTCA], GMM_s.mean(j,:)', mu, options);
                CP   = STMP * GMM_p.cov(:,:,i) * STMP';
                CS   = STMS * GMM_s.cov(:,:,j) * STMS';
            %---------------------------- MCI -----------------------------
            case 'mci'
                % Find TCA (Refined TCA) relative to nominal TCA t0 (match CR3BP reference)
                % GMM mean indexing can be [N x 6]; ensure we pass column vectors
                xp0 = GMM_p.mean(i,:).';
                xs0 = GMM_s.mean(j,:).';
                tTCA = RefineTCA_MCI(xp0, xs0, t0);
                % Propagation only to refined TCA (seconds) – no dense grid needed
                if tTCA ~= t0
                    [~, xP_m] = ode113(ode_fun, [t0, tTCA], xp0, options);
                    [~, xS_m] = ode113(ode_fun, [t0, tTCA], xs0, options);
                    % Propagate covariance via STM (initial->final mapping)
                    [~, STMP] = stateTrans_perturbed_2bp(xp0, tTCA - t0, params, options, t0);
                    [~, STMS] = stateTrans_perturbed_2bp(xs0, tTCA - t0, params, options, t0);
                    % STM maps initial->final: C_f = Phi * C_i * Phi'
                    CP = STMP * GMM_p.cov(:,:,i) * STMP';
                    CS = STMS * GMM_s.cov(:,:,j) * STMS';
                    % enforce symmetry to reduce numerical noise
                    CP = (CP + CP')/2;
                    CS = (CS + CS')/2;
                else %Equal TCA -> New is set equal to old
                    xP_m = xp0.'; % row
                    xS_m = xs0.'; % row
                    CP   = GMM_p.cov(:,:,i);
                    CS   = GMM_s.cov(:,:,j);
                end
                % Convert states from meters to kilometers for Chan_POC
                xP = xP_m(end,:); xS = xS_m(end,:);
                xP(1:3) = xP(1:3) / 1000; xP(4:6) = xP(4:6) / 1000;
                xS(1:3) = xS(1:3) / 1000; xS(4:6) = xS(4:6) / 1000;
                % Convert covariance to km-based units
                CP = CP / 1e6; CS = CS / 1e6;
            otherwise
                error('Unknown AstroMod: %s', AstroMod);
        end

        %Compute GMM POC (match reference behavior: pass full 6x6, Chan extracts 3x3)
        Csum = (CP + CS);
        % enforce symmetry for stability without changing logic
        Csum = (Csum + Csum')/2;
        POC_GMM = POC_GMM + wP(i) * wS(j) * Chan_POC(xP, xS, Csum, Robjs, ChanOrd);
    end
end

end
