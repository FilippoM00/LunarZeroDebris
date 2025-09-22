function cfg = LLO_defaults()
% LLO_DEFAULTS  Centralized default parameter ranges & options for LLO generation/validation.
%
% All numeric ranges are km (altitudes) or degrees (angles) unless specified.
% Modify ONLY here; generation and validation scripts should not hard-code values.
%
% Core fields:
%   alt_range                - Generic altitude band for circular / polar base (km)
%   circular_inc_range       - Inclination sampling for circular family (deg)
%   circular_omega_range     - Argument of periapsis sampling range for circular (deg)
%   circular_RAAN_range      - RAAN sampling range for circular (deg)
%   circular_ecc             - Fixed eccentricity for circular family
%
% Eccentric family:
%   ecc_peri_alt_range       - Periapsis altitude range (km)
%   ecc_apo_alt_range        - Apoapsis altitude range (km)
%   eccentric_inc_range      - Inclination sampling for eccentric family (deg)
%   eccentric_omega_range    - Argument of periapsis sampling range for eccentric (deg)
%   eccentric_RAAN_range     - RAAN sampling range for eccentric (deg)
%
% Polar family:
%   polar_inc_range          - Inclination sampling for polar family (deg)
%   polar_omega_options      - Allowed ω values (deg) (north/south pole alignment)
%   polar_RAAN_range         - RAAN sampling range for polar (deg)
%   polar_circular_probability - Probability a polar sample is circular (0..1)
%   polar_apo_offset_km      - [min max] offset added to peri altitude for eccentric polar apoapsis (km)
%
% Frozen family:
%   frozen_inc_range         - Inclination sampling for frozen family (deg)
%   frozen_omega_range       - Argument of periapsis sampling range for frozen (deg)
%   frozen_RAAN_range        - RAAN sampling range for frozen (deg)
%   frozen_ecc_range         - Eccentricity sampling range
%   frozen_type_weights      - Weights for frozen types [kozai polar nearpolar lowinc highinc] (sum=1)
%
% General switches:
%   seed                     - RNG seed (set [] to skip)
%
cfg = struct();

% Core circular
cfg.alt_range              = [50 400];
cfg.circular_inc_range     = [0 180];
cfg.circular_omega_range   = [0 360];
cfg.circular_RAAN_range    = [0 360];
cfg.circular_ecc           = 0.0;

% Eccentric family - MAXIMUM eccentricity variety
cfg.ecc_peri_alt_range     = [50 200];   % EVEN lower periapsis for higher eccentricity
cfg.ecc_apo_alt_range      = [600 1000]; % EVEN higher apoapsis for maximum eccentricity
cfg.eccentric_inc_range    = [0 180];
cfg.eccentric_omega_range  = [0 360];
cfg.eccentric_RAAN_range   = [0 360];

% Polar specifics - MAXIMUM eccentricity variety
cfg.polar_inc_range        = [85 95];
cfg.polar_omega_options    = [0 180];
cfg.polar_RAAN_range       = [0 360];
cfg.polar_circular_probability = 0.2;    % FURTHER reduced to allow even more eccentric orbits
cfg.polar_apo_offset_km    = [200 1000];  % GREATLY increased from [100 800] for much higher eccentricity

% Frozen sampling - Literature-based families
% Based on Folta & Quinn (2006), Elipe & Lara (2003), Nie & Gurfil (2018)

% Define frozen families explicitly with precise literature values and enhanced parameters
frozen_families.low_inc1 = struct(...
    'inc_deg', 27.05, ...
    'inc_tolerance_deg', 2.0, ...  % ±2° stability boundary
    'ecc_range', [0.0 0.3], ...  % INCREASED from [0.01 0.05] for more variety
    'alt_range_km', [50 200], ...  % Optimal altitude range for stability
    'omega_libration_deg', 5.0, ... % Libration amplitude around ±90°
    'description', 'Low-inclination frozen (J2-Earth balance)', ...
    'perturbation_balance', 'J2_dominant', ... % Primary perturbation mechanism
    'stability_notes', 'Requires careful third-body perturbation balance'...
);

frozen_families.low_inc2 = struct(...
    'inc_deg', 39.23, ...
    'inc_tolerance_deg', 2.0, ...
    'ecc_range', [0.0 0.3], ...  % INCREASED from [0.01 0.05] for more variety
    'alt_range_km', [50 200], ...
    'omega_libration_deg', 5.0, ...
    'description', 'Low-inclination frozen (J2-Earth balance)', ...
    'perturbation_balance', 'J2_dominant', ...
    'stability_notes', 'Requires careful third-body perturbation balance'...
);

frozen_families.magic1 = struct(...
    'inc_deg', 69.00, ...
    'inc_tolerance_deg', 3.0, ...  % ±3° stability boundary
    'ecc_range', [0.00 0.3], ...
    'alt_range_km', [50 500], ...  % Broader stability range
    'omega_libration_deg', 3.0, ... % Smaller libration due to natural stability
    'description', 'Magic inclination (natural stability)', ...
    'perturbation_balance', 'natural_stability', ...
    'stability_notes', 'Naturally stable against multiple perturbations'...
);

frozen_families.magic2 = struct(...
    'inc_deg', 110.96, ...
    'inc_tolerance_deg', 3.0, ...
    'ecc_range', [0.00 0.3], ...
    'alt_range_km', [50 500], ...
    'omega_libration_deg', 3.0, ...
    'description', 'Magic inclination (natural stability)', ...
    'perturbation_balance', 'natural_stability', ...
    'stability_notes', 'Naturally stable against multiple perturbations'...
);

frozen_families.polar = struct(...
    'inc_deg', 90.00, ...
    'inc_tolerance_deg', 1.0, ...  % ±1° stability boundary
    'ecc_range', [0.0 0.2], ...  % INCREASED from [0.01 0.08] for more variety
    'alt_range_km', [50 300], ...  % Critical eccentricity varies with altitude
    'omega_libration_deg', 5.0, ...
    'description', 'Polar frozen (conservative eccentricity)', ...
    'perturbation_balance', 'J2_conservative', ...
    'stability_notes', 'Uses conservative eccentricity range for stability'...
);

cfg.frozen_families = frozen_families;

% Frozen family weights (should sum to 1)
cfg.frozen_family_weights = [0.2, 0.2, 0.2, 0.2, 0.2]; % Equal weighting for all families

% Legacy frozen parameters (for backward compatibility) - UPDATED ranges
cfg.frozen_inc_range       = [0 180];
cfg.frozen_omega_range     = [0 360];
cfg.frozen_RAAN_range      = [0 360];
cfg.frozen_ecc_range       = [0.001 0.12];  % UPDATED from [0.001 0.05] to match new ranges
cfg.frozen_type_weights    = [0.2 0.2 0.2 0.2 0.2]; % Order: kozai, polarfrozen, nearpolar, lowinc, highinc

% General LLO family - HIGH eccentricity for stability analysis
cfg.general_inc_range      = [10 170];      % Broad inclination range
cfg.general_omega_range    = [0 360];
cfg.general_RAAN_range     = [0 360];
cfg.general_ecc_range      = [0.1 0.4];     % INCREASED from [0.05 0.3] for higher eccentricity
cfg.general_alt_range      = [50 400];      % Reasonable altitude range

% HIGHLY Eccentric LLO family - For extreme cases
cfg.high_ecc_inc_range     = [0 180];       % Full inclination range
cfg.high_ecc_omega_range   = [0 360];
cfg.high_ecc_RAAN_range    = [0 360];
cfg.high_ecc_peri_range    = [50 50];       % Very low periapsis
cfg.high_ecc_apo_range     = [800 2000];    % Very high apoapsis
cfg.high_ecc_alt_range     = [50 100];      % Use periapsis as reference

% --- PHYSICAL AND CR3BP CONSTANTS ---
% These should not be modified unless you are changing the fundamental model.

% Primary gravitational parameters (m^3/s^2)
params.mu_earth = 3.986004418e14;
params.GM_moon = 4.902800122e12;  % Updated to GRAIL GRGM1200A precision
params.GM_earth = params.mu_earth;

% Moon's mean radius (meters)
params.R_moon = 1737.4 * 1000;

% Moon's mean radius (km)
params.R_moon_km = params.R_moon / 1000;

% Earth-Moon system parameters
params.D_em = 384400 * 1000;        % Mean Earth-Moon distance (m)
params.w_moon = 2.6617e-6;          % Mean rotation rate of the Moon (rad/s)
params.n_earth = 2.6617e-6;         % Mean motion of Earth around Moon (rad/s) - same as lunar rotation

% Lunar gravitational harmonics (un-normalized, from GRAIL GRGM1200A model)
params.J2 = 2.0339e-4;   % Updated to GRGM1200A
params.J3 = -9.2e-6;     % Updated to GRGM1200A
params.C22 = 2.243e-5;   % Updated to GRGM1200A

% CR3BP characteristic quantities
L_star = 384400 * 1000; % Characteristic length (Earth-Moon distance in m)
T_star = 27.321661 * 86400; % Characteristic time (sidereal month in s)
V_star = L_star / (T_star / (2*pi)); % Characteristic velocity

params_cr3bp.mu = params.GM_moon / (params.mu_earth + params.GM_moon);
params_cr3bp.L_star = L_star;
params_cr3bp.V_star = V_star;
params_cr3bp.T_star = T_star;
params_cr3bp.R_moon_n = params.R_moon / L_star; % Moon radius, non-dimensional

% Add the structs to the main cfg output
cfg.params = params;
cfg.params_cr3bp = params_cr3bp;

end
