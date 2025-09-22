function [C, Ceci] = CovGen(err, seed, k)
%
% Object Covariance Matrix Generator
%
%DESCRIPTION:
%This code generates the covariance matrix of an orbiting object.
%
%PROTOTYPE
%   [C] = CovGen(err)
%   [C] = CovGen(err, seed)
%   [C] = CovGen(err, seed, k)
%   [C, Ceci] = CovGen(err)
%   [C, Ceci] = CovGen(err, seed)
%   [C, Ceci] = CovGen(err, seed, k)
%
%--------------------------------------------------------------------------
% INPUTS:
%   cov        [1x1]       Mean Covariance             [km]
%   seed       [1x1]       RNG Seed                    [-] (optional)
%   k          [1x1]       Scaling Factor              [-] (optional)
%--------------------------------------------------------------------------
% OUTPUTS:
%   C          [nxn]       Covariance Matrix           [-] (see NOTES)
%   Ceci       [nxn]       Covariance Matrix (BCI)     [-] (see NOTES)
%--------------------------------------------------------------------------
%
%NOTES:
% - The input "err" can be:
%     [1x1] = Position or Velocity  -> [3x3] Cov. Mat.
%     [1x2] = Position and Velocity -> [6x6] Cov. Mat.
%   Note also that the uncertainties in "err" now shall be not-squared!
%
% - The "seed" input shall be an integer number, which will set the RNG
%   working condition.
%
% - The Scaing Factor "k" can be given as input to correct the covariance
%   in order to increase or decrease the covariance.
%
% - To verify the correctness of the generation, check that eigenvalues of
%   C are all >=0 and in the order of tens of km. Moreover, recall that in
%   the RTM frame, the T-semimajor axis shall be the longest.
%
% - The outputs "C" and "Ceci" have to be intended in this way:
%   * if the input is meant for Synodic use (errors are dimensionless) then C
%     is Synodic and the (optional) output Ceci is the
%     re-dimensionalization in EBCI of C.
%   * if the input is meant for BCI use (errors are dimensional) then C is
%     directly in BCI and Ceci shall not be taken.
%
%CALLED FUNCTIONS:
% (none)
%
%UPDATES:
% 2022/06/21, Luigi De Maria - Expanded the function to position and
%             velocity case, changed from cov^2 to cov in input
% 2022/09/13, Luigi De Maria - Added the ECi version of C in Synodic case
%
%REFERENCES:
% (none)
%
%AUTHOR(s):
%Luigi De Maria, 2022
%

%% Main Code

%Random Number Generator Seed Definition
if nargin > 1
	rng(seed);
end

%Random Matrix Generation
if     length(err) == 1
	C = rand(3);
elseif length(err) == 2
	C = rand(6);
end

%Make it Semi-Definite Positive
C = C * C';

%Scaled by Covariance
if     length(err) == 1
	C = C .* err^2;
elseif length(err) == 2
	C(1:3,1:3) = C(1:3,1:3) .* err(1)^2;
	C(1:3,4:6) = C(1:3,4:6) .* (err(1)*err(2));
	C(4:6,1:3) = C(4:6,1:3) .* (err(2)*err(1));
	C(4:6,4:6) = C(4:6,4:6) .* err(2)^2;
end

%Scale by Scaling Factor
if nargin == 3
	C = C .* k;
end

%Convertion to BCI (re-dimensionalization)
if nargout == 2
	%Synodic Constants
	tau = 4.34811305;               %Synodic Time Constant [days]
	D   = 384405;                   %Synodic Characteristic Distance [km]
	%Re-dimensionalization
	if     length(err) == 1
		Ceci = C .* D^2;
	elseif length(err) == 2
		Ceci(1:3,1:3) = C(1:3,1:3) .* D^2;
		Ceci(1:3,4:6) = C(1:3,4:6) .* D*(D/(tau*3600*24));
		Ceci(4:6,1:3) = C(4:6,1:3) .* D*(D/(tau*3600*24));
		Ceci(4:6,4:6) = C(4:6,4:6) .* (D/(tau*3600*24))^2;
	end
end

%Reset RNG
rng('default');

end
