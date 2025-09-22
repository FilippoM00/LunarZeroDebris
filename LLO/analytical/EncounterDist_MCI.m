function [min_dist, tEncounter] = EncounterDist_MCI(state, tspan, params, Method, Tsearch)
%
% Encounter Minimum Distance Function
%
%DESCRIPTION:
%This code provides the function for the computation of the encounter
%minimum distance.
%
%PROTOTYPE
%   [min_dist] = EncounterDist_MCI(state, tspan, params)
%   [min_dist] = EncounterDist_MCI(state, tspan, params, Method)
%   [min_dist] = EncounterDist_MCI(state, tspan, params, Method, Tsearch)
%   [min_dist, tEncounter] = EncounterDist_MCI(state, tspan, params)
%   [min_dist, tEncounter] = EncounterDist_MCI(state, tspan, params, Method)
%   [min_dist, tEncounter] = EncounterDist_MCI(state, tspan, params, Method, Tsearch)
%
%--------------------------------------------------------------------------
% INPUTS:
%   state      [1x12]      State Vector              [m],[m][s-1]
%   tspan      [1xn]       Time Span Vector          [s]
%   params     [struct]    Physical parameters for MCI perturbed 2BP [-]
%   Method     [1x1]       Interpolation Method      [-] (string) (NOTES)
%   Tsearch    [1x2]       Search Time Range         [s] (optional)
%--------------------------------------------------------------------------
% OUTPUTS:
%   dmin_distx [1x1]       Minimum Distance          [m]
%   tEncounter [1x1]       Encounter Time Instant    [s] (or [-]) (optional)
%--------------------------------------------------------------------------
%
%NOTES:
% - The State Vector shall contain both Primary and Secondary states.
% - This MCI version mirrors the CR3BP pipeline but always uses the
%   perturbed 2-body dynamics in the Moon-Centered Inertial frame.
% - The input "Method" defines the interpolation method adopted.
%   * Method = 'polyfit3' --> Polynomial fitting 3rd Order
%   * Method = 'spline'   --> Interpolation w/splines (NOT ADVICED)
%   * Method = 'phillips' --> Linear interpolation of [1] (NOT ADVICED)
%   * Method = 'refine'   --> Propagation Refinement
% - When using the 'refine' method be careful of the time step used in
%   "tspan". The refinement considers 2 points and creates a linspace time
%   vector of 300 points (default value).
%   Therefore, it is adviced to have:
%     - Maximum Time Step of tspan = 10s
%     - Time Step of Refined tspan < 0.1s
%   As a tip: consider the velocities involved and the time step that you
%   would have in the refinement. If v*dt gives a dx in the order of cm
%   it's ok, otherwise you may have errors of meters and
%   overestimate/underestimate the number of collisions.
%
%CALLED FUNCTIONS:
% two_body_perturbed_rhs
%
%UPDATES:
% 2025/09/15, EncounterDist_MCI rewritten to mirror CR3BP pipeline with MCI dynamics.
%
%REFERENCES:
% [1] "Spacecraft Collision Probability Estimation for Rendezvous and
%      Proximity Operations". Michael R. Phillips.
%
%AUTHOR(s):
%Luigi De Maria (original idea), MCI adaptation 2025
%

%% Main Code

%Constants and Parameters
ode_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13); %ODE Settings
ode_fun     = @(t,s) two_body_perturbed_rhs(t, s, params);

% Provide default method
if nargin < 4 || isempty(Method)
    Method = 'refine';
end
% Provide default search window if not supplied
if nargin < 5 || isempty(Tsearch)
    Tsearch = [tspan(1), tspan(end)];
end

%Split State Vectors
xP = state(1:6);
xS = state(7:12);

%------------------------- Propagator (MCI) ------------------------------
%Extract Subinterval Indices (inclusive bounds)
indx  = find(tspan > (Tsearch(1)) & tspan < (Tsearch(2)));
%Refinement Time Step [s] (set by dtref variable; value is in seconds)
dtref = 1e-3;
%Propagation
if isrow(xP) && isrow(xS)
    [~, state_prop_P] = ode45(ode_fun, tspan, xP', ode_options);
    [~, state_prop_S] = ode45(ode_fun, tspan, xS', ode_options);
elseif iscolumn(xP) && iscolumn(xS)
    [~, state_prop_P] = ode113(ode_fun, tspan, xP,  ode_options);
    [~, state_prop_S] = ode113(ode_fun, tspan, xS,  ode_options);
else
    error('Make sure ICs are both row or column vectors!')
end

%----------------------- Minimum Distance Search --------------------------
%Reduced States (in Search Window)
stateP_red = state_prop_P(indx,1:6);
stateS_red = state_prop_S(indx,1:6);

%Relative Position and Velocity
Dr = stateP_red(:,1:3)-stateS_red(:,1:3);
Dv = stateP_red(:,4:6)-stateS_red(:,4:6);

%DrDv Condition
DrDv = dot(Dr,Dv,2);

%New Sub-Time Span
tspan_new = tspan(indx);

%Zero-Crossings Detection
aux = diff(sign(DrDv),1,1);
%Zero-Crossings Positions
zeroCrossings = find(aux);

%Zero-Crossings Detected
if ~isempty(zeroCrossings)
    %Dummy Value of Relative Distance
    RelDistInterp = inf;
    for k = 1 : length(zeroCrossings)
        %---------------------- Polynomial Fitting ------------------------
        if     strcmp(Method,'polyfit3')
            %Adjust Interpolation Samples if Sign Change Location near Boundaries
            if     (zeroCrossings(k) == 1)
                lbound = zeroCrossings(k);
                rbound = zeroCrossings(k) + 3;
            elseif (zeroCrossings(k) == 2)
                lbound = zeroCrossings(k) - 1;
                rbound = zeroCrossings(k) + 2;
            elseif (zeroCrossings(k) == length(indx)-1)
                lbound = zeroCrossings(k) - 2;
                rbound = zeroCrossings(k) + 1;
            elseif (zeroCrossings(k) == length(indx))
                lbound = zeroCrossings(k) - 3;
                rbound = zeroCrossings(k);
            else
                lbound = zeroCrossings(k) - 2;
                rbound = zeroCrossings(k) + 1;
            end

            %Polynomial Fitting of 3rd Order (DrDv)
            [p1,~,~] = polyfit(tspan_new(lbound:rbound)', DrDv(lbound:rbound,:), 3);
            %Polynomial Fitting of 3rd Order (Dr)
            [p2,~,mu2] = polyfit(tspan_new(lbound:rbound)', vecnorm(Dr(lbound:rbound,1:3),2,2), 3);
            %Roots of DrDv Polynomial
            r  = roots(p1);
            r  = real(r(imag(r)==0)); %Save only real roots
            
            %Determine if Roots are in Time Bounds
            %aux = find((r > tspan_new(lbound)) & (r < tspan_new(rbound)));
            for kk = 1 : length(r)
               %Find Relative Distances w/Roots of DrDv
                Dreval = polyval(p2,r(kk),[],mu2);
                %Collision Condition (Minimum of Evaluated Dr)
                if Dreval < RelDistInterp
                    RelDistInterp = Dreval;
                    %Save Encounter Time Instant
                    if nargout == 2
                        tEncounter = r(kk);
                    end
                end
            end
            
        %------------------------ Spline Fitting --------------------------
        elseif strcmp(Method,'spline')
            %Adjust Interpolation Samples if Sign Change Location near Boundaries
            if     (zeroCrossings(k) == 1)
                lbound = zeroCrossings(k);
                rbound = zeroCrossings(k) + 2;
            elseif (zeroCrossings(k) == length(indx))
                lbound = zeroCrossings(k) - 2;
                rbound = zeroCrossings(k);
            else
                lbound = zeroCrossings(k) - 1;
                rbound = zeroCrossings(k) + 1;
            end

            %Polynomial Fitting of 3rd Order (Dr)
            Drnorm = vecnorm(Dr(lbound:rbound,1:3),2,2);

            %Spline Interpolation
            pp = spline(tspan_new(lbound:rbound)',DrDv(lbound:rbound,:));
            zeros_Spline = fnzeros(pp);

            for kk = 1 : length(zeros_Spline)
               %Find Relative Distances w/Roots of DrDv
                Dreval = spline(tspan_new(lbound:rbound)',Drnorm,zeros_Spline(kk));
                %Collision Condition (Minimum of Evaluated Dr)
                if Dreval < RelDistInterp
                    RelDistInterp = Dreval;
                    %Save Encounter Time Instant
                    if nargout == 2
                        tEncounter = zeros_Spline(kk);
                    end
                end
            end
            
        %----------------------- External Function ------------------------
    elseif strcmp(Method,'phillips')
            %Identity Matrix
            I3 = eye(3);
            for j = 1 : length(tspan_new)-1
                %Relative Pos. @t_i (6 elements, we only need position part but keep as vector)
                ri = stateP_red(j,:) - stateS_red(j,:);
                %Relative Pos. @t_{i+1}
                ri_next = stateP_red(j+1,:) - stateS_red(j+1,:);
                %Vector from i to i+1 Rel.Pos.
                i_t = (ri - ri_next) / norm(ri - ri_next);
                %Linearly Interpolated Relative Pos. (Mid-Point)
                rn = (I3 - i_t'*i_t) * ri';
                %Norm of rn
                Dreval = norm(rn);
                %Check Minimum Rel.Dist.
                if Dreval < RelDistInterp
                    RelDistInterp = Dreval;
                    %Save Encounter Time Instant
                    if nargout == 2
                        tEncounter = tspan_new(j);
                    end
                end
            end
            
        %-------------------- Propagation Refinement ----------------------
        elseif strcmp(Method,'refine')
            %Reduce Boundaries of Refined Propagation
            if   (zeroCrossings(k) == length(indx))
                lbound = zeroCrossings(k) - 1;
                rbound = zeroCrossings(k);
            else
                lbound = zeroCrossings(k);
                rbound = zeroCrossings(k) + 1;
            end
            %Refined Time Span
            if     tspan_new(lbound) < tspan_new(rbound)
                tspan_ref = [tspan_new(lbound) :  dtref : tspan_new(rbound)];
            elseif tspan_new(lbound) > tspan_new(rbound)
                tspan_ref = [tspan_new(lbound) : -dtref : tspan_new(rbound)];
            end
            %Propagate Objects in Refine Time Span
            [~, stateP_ref] = ode113(ode_fun, tspan_ref, stateP_red(lbound,:)', ode_options);
            [~, stateS_ref] = ode113(ode_fun, tspan_ref, stateS_red(lbound,:)', ode_options);
            %Norm of Relative Distance in Refined Time Span
            Dreval = vecnorm(stateP_ref(:,1:3)-stateS_ref(:,1:3),2,2);
            %Check Minimum Distance Improvement
            [min_val_ref, min_idx_ref] = min(Dreval);
            if min_val_ref < RelDistInterp
                %Minimum Relative Distance
                RelDistInterp = min_val_ref;
                %Save Encounter Time Instant (scalar)
                if nargout == 2
                    tEncounter = tspan_ref(min_idx_ref);
                end
            end
            
        end
    end
else
    [RelDistInterp, min_idx_coarse] = min(vecnorm(Dr,2,2));
    if nargout == 2
        tEncounter = tspan_new(min_idx_coarse);
    end
end

%Minimum of MD
min_dist = RelDistInterp;
% Safety clamp: ensure tEncounter lies within the requested window
% (No additional clamping in reference implementation)

end
