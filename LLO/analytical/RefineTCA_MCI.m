function [TCAnew, ckTCA] = RefineTCA_MCI(xP, xS, TCAold, Tback, Tfor)
%
% TCA Refinement for MCI Perturbed 2-Body Dynamics
%
%DESCRIPTION:
%This code provides the implementation of TCA refinement routine for MCI
%perturbed two-body dynamics (Moon gravity + Earth 3rd-body + J2).
%
%Adapted from RefineTCA.m (CR3BP) to MCI framework.
%
%PROTOTYPE
%   [TCAnew] = RefineTCA_MCI(xP, xS, TCAold)
%   [TCAnew] = RefineTCA_MCI(xP, xS, TCAold, Tback, Tfor)
%   [TCAnew, ckTCA] = RefineTCA_MCI(xP, xS, TCAold)
%   [TCAnew, ckTCA] = RefineTCA_MCI(xP, xS, TCAold, Tback, Tfor)
%
%--------------------------------------------------------------------------
% INPUTS:
%   xP         [6x1]       Primary State @TCAold           [m, m/s]
%   xS         [6x1]       Secondary State @TCAold         [m, m/s]
%   TCAold     [1x1]       Old TCA                         [s]
%   Tback      [1x1]       Left-end Research Time Window   [s] (optional)
%   Tfor       [1x1]       Right-end Research Time Window  [s] (optional)
%--------------------------------------------------------------------------
% OUTPUTS:
%   TCAnew     [1x1]       Refined TCA                     [s]
%   ckTCA      [1x1]       DrDv check at TCA               [m^2/s] (optional)
%--------------------------------------------------------------------------
%
%NOTES:
% - Both TCAold and TCAnew are in "seconds".
% - The research time span is fixed to [-2hours,+2hours] around TCAold.
%   Use optional inputs if another window is desired.
% - This function also plots DrDv and Dr trends with time (commented out).
% - Use "ckTCA" for output on command window of TCA data.
%
%AUTHOR(s):
% Adapted by GitHub Copilot, 2025
% Original by Luigi De Maria, 2022
%

%% Main Code

%------------------------------ Quantities --------------------------------

% --- Input validation ---
if nargin ~= 3 && nargin ~= 5
    error('RefineTCA_MCI:InvalidInputs', 'Expected 3 or 5 inputs: (xP, xS, TCAold [, Tback, Tfor])');
end
if ~isvector(xP) || ~isvector(xS) || numel(xP) ~= 6 || numel(xS) ~= 6
    error('RefineTCA_MCI:StateSize', 'xP and xS must be 6-element vectors');
end
% Ensure column vectors
xP = xP(:);
xS = xS(:);
if ~isscalar(TCAold) || ~isfinite(TCAold)
    error('RefineTCA_MCI:TCAold', 'TCAold must be a finite scalar in seconds');
end
if nargin == 5
    if ~isscalar(Tback) || ~isscalar(Tfor) || Tback < 0 || Tfor < 0
        error('RefineTCA_MCI:Windows', 'Tback and Tfor must be nonnegative scalars in seconds');
    end
end

% Load parameters (assume LLO_defaults available)
cfg = LLO_defaults();
params = cfg.params;

% ODE settings
ode_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
ode_fun = @(t,x) two_body_perturbed_rhs(t, x, params); % MCI ODE (replaces CR3BP ODE)

% Number of Points
n = 30000;

% Time Spans (in seconds)
if nargin == 3 % Default
    tspanF = linspace(TCAold, TCAold + 2*3600, floor(n/2)); % Forward [s]
    tspanB = linspace(TCAold, TCAold - 2*3600, floor(n/2)); % Backward [s]
elseif nargin == 5 % User-defined
    tspanF = linspace(TCAold, TCAold + Tfor, floor(n/2)); % Forward [s]
    tspanB = linspace(TCAold, TCAold - Tback, floor(n/2)); % Backward [s]
else
    error('Wrong number of inputs.')
end

% Propagation (MCI)
[~, xP_F] = ode113(ode_fun, tspanF, xP, ode_options); % Primary Forward
[~, xP_B] = ode113(ode_fun, tspanB, xP, ode_options); % Primary Backward
[~, xS_F] = ode113(ode_fun, tspanF, xS, ode_options); % Secondary Forward
[~, xS_B] = ode113(ode_fun, tspanB, xS, ode_options); % Secondary Backward

% Flip backward arrays
xP_B = flip(xP_B, 1);
xS_B = flip(xS_B, 1);
tspanB = flip(tspanB, 2);

% Concatenate [-t, t]
x1 = [xP_B(1:end-1,:); xP_F];
x2 = [xS_B(1:end-1,:); xS_F];
tsim = [tspanB(1:end-1), tspanF];

%{
% PLOTTING (commented out)
ck = dot((x1(:,1:3)-x2(:,1:3)), (x1(:,4:6)-x2(:,4:6)), 2);
reldist = vecnorm(x1(:,1:3)-x2(:,1:3), 2, 2);

% DrDv plot
figure; semilogy(tsim, abs(ck), 'b', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('DrDv [m^2/s]'); title('DrDv Condition'); grid on;

% Relative Distance plot
figure; semilogy(tsim, reldist, 'b', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Relative Distance [m]'); title('Relative Distance'); grid on;
%}

%---------------------------- TCA Refinement ------------------------------

% Relative Position and Velocity
Dr = x1(:,1:3) - x2(:,1:3);
Dv = x1(:,4:6) - x2(:,4:6);

% DrDv Condition
DrDv = dot(Dr, Dv, 2);

% Zero-Crossings Detection
aux = diff(sign(DrDv), 1, 1);
zeroCrossings = find(aux);

% Zero-Crossings Detected
if ~isempty(zeroCrossings)
    % Dummy value
    RelDistInterp = inf;
    for k = 1:length(zeroCrossings)
        % Reduce boundaries
        if zeroCrossings(k) == length(tsim)
            lbound = zeroCrossings(k) - 1;
            rbound = zeroCrossings(k);
        else
            lbound = zeroCrossings(k);
            rbound = zeroCrossings(k) + 1;
        end
        % Refined time span
        if tsim(lbound) < tsim(rbound)
            tspan_ref = [tsim(lbound) :  1e-5 : tsim(rbound)]; % 0.001 s steps (match CR3BP effective seconds step)
        elseif tsim(lbound) > tsim(rbound)
            tspan_ref = [tsim(lbound) : -1e-5 : tsim(rbound)];
        end
        % Propagation with refined span
        [~, stateP_ref] = ode113(ode_fun, tspan_ref, x1(lbound,:)', ode_options);
        [~, stateS_ref] = ode113(ode_fun, tspan_ref, x2(lbound,:)', ode_options);
        % Relative distance
        Dreval = vecnorm(stateP_ref(:,1:3) - stateS_ref(:,1:3), 2, 2);
        % Check min distance
        if min(Dreval) < RelDistInterp
            RelDistInterp = min(Dreval);
            cand = tspan_ref(Dreval == min(Dreval));
            % If multiple equal minima, pick the one closest to TCAold
            if numel(cand) > 1
                [~, idxmin] = min(abs(cand - TCAold));
                TCAnew = cand(idxmin);
            else
                TCAnew = cand;
            end
        end
    end
else
    % Min distance
    [RelDistInterp, minIndx] = min(vecnorm(Dr, 2, 2));
    TCAnew = tsim(minIndx);
end

% Ensure scalar TCAnew (CR3BP reference uses scalar). If multiple remain, pick closest to TCAold.
if numel(TCAnew) > 1
    [~, idxmin] = min(abs(TCAnew - TCAold));
    TCAnew = TCAnew(idxmin);
end

%Check Final Time
if TCAnew == TCAold
    TCAnew   = TCAold + eps;
end

% TCA Data
if nargout == 2
    % Propagation to refined TCA
    [~, x_TCArefP] = ode113(ode_fun, [TCAold, TCAnew], xP, ode_options);
    [~, x_TCArefS] = ode113(ode_fun, [TCAold, TCAnew], xS, ode_options);
    
    % DrDv check
    ckTCA = dot((x_TCArefP(end,1:3) - x_TCArefS(end,1:3)), (x_TCArefP(end,4:6) - x_TCArefS(end,4:6)));
    
    % Print results (match CR3BP messages)
    disp('The refinement process determined the cloasest approach at: [s]')
    disp(TCAnew)
    if (TCAnew == tsim(1)) || (TCAnew == tsim(end))
        warning('The found TCA is near boundaries of tspan, maybe too narrow search window. Consider enlarging it.')
    end
    disp('The check on DrDv at TCA: [-]')
    disp(ckTCA)
    disp('The resulting relative distance at TCA: [m]')
    disp(RelDistInterp)
end

end
