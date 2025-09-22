function [Rb, covariance_bp, relPos_bp] = compute_bplane(mean_state_1,...
    covariance_pos_1,mean_state_2,covariance_pos_2,varargin)
% Compute the main features of Bplane.
% -------------------------------------------------------------------------
% INPUTS:
% - mean_state_1    : primary object mean state in ECI reference frame
% - covariance_pos_1: position covariance of the primary object in ECI
%                     reference frame
% - mean_state_2    : secondary object mean state in ECI reference frame
% - covariance_pos_2: position covariance of the secondary object in ECI
%                     reference frame
% - Kp              : primary object covariance sensitivity parameter
%                     (optional)
% - Ks              : secondary object covariance sensitivity parameter
%                     (optional)
% -------------------------------------------------------------------------
% OUTPUT:
% - Rb           : rotation matrix from ECI to Bplane
% - covariance_bp: global covariance in the Bplane reference frame
% - relPos_bp    : relative position between primary and secondary in the
%                  Bplane reference frame
% -------------------------------------------------------------------------
% Author:  Marco Felice Montaruli, Politecnico di Milano, 02 February 2022
%          e-mail: marcofelice.montaruli@polimi.it

% Minimum input number
n_input_min = 4;

% Default Kp-Ks values
Kp = 1;
Ks = 1;


% Record input covariance matrix and print_tag if given as optional input
if nargin>n_input_min
    for index=1:(nargin-n_input_min)
        if ~ischar(varargin{index}) && ~isstring(varargin{index})
            continue
        end
        switch varargin{index}
            case 'Kp'
                Kp = varargin{index+1};
            case 'Ks'
                Ks = varargin{index+1};
            otherwise
                ME = MException('PoliMiSST:ConjunctionAnalysis:invalidOptionalInput', ...
                    'The optional input is not valid');
                throw(ME)
        end
    end
end


% Matrix to rotate in Bplane reference frame
Rb = R_eci2bplane(mean_state_1(4:6)',mean_state_2(4:6)');

% Objects mean position
r_1_eci = mean_state_1(1:3); % [km]
r_2_eci = mean_state_2(1:3); % [km]

% Mean relative distance between the two objects in ECI
relPos_eci = r_2_eci - r_1_eci; % [km] ECI r.f.

% Mean relative distance between the two objects in bplane
relPos_bp_3D = Rb*relPos_eci'; % [km]  B-Plane r.f.
relPos_bp = [relPos_bp_3D(1); 0; relPos_bp_3D(3)];


% Global covariance projected in b-plane
covariance_bp_1 = Rb*covariance_pos_1*Rb'; %[km^2]
covariance_bp_2 = Rb*covariance_pos_2*Rb'; %[km^2]
covariance_bp   = Kp*covariance_bp_1+Ks*covariance_bp_2; %[km^2]

% Covariance 2D
covariance_bp = [covariance_bp(1,1), covariance_bp(1,3); ...
    covariance_bp(3,1), covariance_bp(3,3)]; %[km^2]


end
