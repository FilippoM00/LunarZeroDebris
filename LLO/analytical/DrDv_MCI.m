function [DRDV] = DrDv_MCI(t0, t, xP, xS, options, AstroMod, scale)
%
% DeltaR-DeltaV Function
%
%DESCRIPTION:
%This code provides the function for the computation of dot-product between
%DeltaR and DeltaV. It is meant to be used as handle-function for the TCA
%refinement process.
%
%PROTOTYPE
%   [DRDV] = DrDv(t0, t, xP, xS, options, AstroMod, scale)
%
%--------------------------------------------------------------------------
% INPUTS:
%   t0         [1x1]       Initial Time          [s]
%   t          [1x1]       Simulation Time       [s]
%   xP         [6x1]       Primary IC            [km,km/s] for 'kep' | [m,m/s] for 'mci'
%   xS         [6x1]       Secondary IC          [km,km/s] for 'kep' | [m,m/s] for 'mci'
%   options    [---]       ODE Options           [-]
%   AstroMod   [1x1]       Astrodynamical Model  [string] ('kep' or 'mci')
%   scale      [1x1]       Scaling Factor        [-]
%--------------------------------------------------------------------------
% OUTPUTS:
%   DRDV       [1x1]       Condition Check       [km^2/s] for 'kep' | [m^2/s] for 'mci'
%--------------------------------------------------------------------------
%
%NOTES:
% - The input "t0" shall be the nominal TCA of the distribution. It is used
%   to refine the TCA.
% - The input "AstroMod" shall be 'kep' or 'mci', on the base of the
%   problem under study.
% - The use of "scale" is adviced when dealing with optimization.
%
%CALLED FUNCTIONS:
% analytic_kep_prop
%
%UPDATES:
% (none)
%
%REFERENCES:
% (none)
%
%AUTHOR(s):
%Luigi De Maria, 2022
%

%% Main Code

%Compute DrDv
switch AstroMod
	%KEPLERIAN
	case 'kep'
		%Propagation
		mu   = 4902.800076; %Moon Planetary Parameter [km3][s-2]
		[x1] = analytic_kep_prop(xP', [t0 t], mu);
		[x2] = analytic_kep_prop(xS', [t0 t], mu);

		%Check Condition
		if nargin == 6
			scale = 1;
		end
		DRDV = dot((x1(1:3)-x2(1:3)), (x1(4:6)-x2(4:6))) * scale;
	%MCI perturbed two-body (replaces CR3BP case)
	case 'mci'
		% Load params from defaults to keep signature identical to CR3BP
		cfg = LLO_defaults();
		params = cfg.params;
		%Propagation (MCI)
		[~, x1] = ode113(@(tt,xx) two_body_perturbed_rhs(tt, xx, params), [t0 t], xP, options);
		[~, x2] = ode113(@(tt,xx) two_body_perturbed_rhs(tt, xx, params), [t0 t], xS, options);

		%Check Condition
		if nargin == 6
			scale = 1;
		end
		DRDV = dot((x1(end,1:3)-x2(end,1:3)), (x1(end,4:6)-x2(end,4:6))) * scale;
	otherwise
		error('DrDv_MCI:InvalidModel','AstroMod must be ''kep'' or ''mci''');
end

end
