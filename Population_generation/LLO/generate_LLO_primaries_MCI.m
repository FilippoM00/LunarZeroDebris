function primaries_table = generate_LLO_primaries_MCI(n_samples, family_type, alt_range, cfg_override)
% GENERATE_LLO_PRIMARIES_MCI - Generate LLO primaries in a Moon-Centered Inertial frame.
%
% This function creates LLO initial conditions for a specific orbit family,
% outputting the state vector in SI units for use with a perturbed two-body
% dynamics model.
%
% For frozen orbits, implements literature-based families:
% - Low-inclination frozen (27°, 39°): J2-Earth perturbation balance
% - Magic inclinations (69°, 111°): Natural stability against perturbations
% - Polar frozen (90°): Critical eccentricity for argument of perigee libration
%
% Inputs:
%   n_samples     - Number of LLO primaries to generate
%   family_type   - 'Circular' | 'Eccentric' | 'Polar' | 'Frozen'
%   alt_range     - [min_alt max_alt] km above lunar surface.
%   cfg_override  - (optional struct) fields to override defaults.
%
% Outputs:
%   primaries_table - Table with columns including state vectors (x,y,z,vx,vy,vz)
%                     in meters and m/s in the MCI frame.
%
% References:
%   Folta & Quinn (2006) - Lunar frozen orbits
%   Elipe & Lara (2003) - Frozen orbits about the Moon
%   Nie & Gurfil (2018) - Lunar frozen orbits revisited

% Load defaults and constants
cfg = LLO_defaults();
if nargin >= 4, cfg = apply_overrides(cfg, cfg_override); end
if nargin < 3 || isempty(alt_range), alt_range = cfg.alt_range; end

params = cfg.params;

% Respect external RNG seeding. If cfg.seed provided via overrides, apply it once here.
if isfield(cfg, 'seed') && ~isempty(cfg.seed) && isreal(cfg.seed)
    rng(cfg.seed);
    %fprintf('Using configured RNG seed for %s family generation (seed=%g)\n', family_type, cfg.seed);
end

% Pre-allocate arrays
IDs = (1:n_samples)';
Families = cell(n_samples, 1);
Altitudes = zeros(n_samples, 1);
Inclinations = zeros(n_samples, 1);
Eccentricities = zeros(n_samples, 1);
ArgPeriapsis = zeros(n_samples, 1);
RAANs = zeros(n_samples, 1);
TrueAnomalies = zeros(n_samples, 1);
X = zeros(n_samples, 1); Y = zeros(n_samples, 1); Z = zeros(n_samples, 1);
VX = zeros(n_samples, 1); VY = zeros(n_samples, 1); VZ = zeros(n_samples, 1);

% Generation loop
for i = 1:n_samples
    % --- This section is largely identical to the CR3BP generator ---
    % --- It samples orbital elements based on the family type ---
    
    % Sample altitude and eccentricity based on family
    switch lower(family_type)
        case 'circular'
            alt_km = alt_range(1) + (alt_range(2) - alt_range(1)) * rand();
            ecc = cfg.circular_ecc;
            inc_deg = cfg.circular_inc_range(1) + diff(cfg.circular_inc_range) * rand();
            omega_deg = cfg.circular_omega_range(1) + diff(cfg.circular_omega_range) * rand();
            Omega_deg = cfg.circular_RAAN_range(1) + diff(cfg.circular_RAAN_range) * rand();
        case 'eccentric'
            peri_alt_km = cfg.ecc_peri_alt_range(1) + diff(cfg.ecc_peri_alt_range) * rand();
            apo_alt_km = cfg.ecc_apo_alt_range(1) + diff(cfg.ecc_apo_alt_range) * rand();
            if apo_alt_km <= peri_alt_km, apo_alt_km = peri_alt_km + 1; end
            r_p = (params.R_moon_km + peri_alt_km) * 1000;
            r_a = (params.R_moon_km + apo_alt_km) * 1000;
            ecc = (r_a - r_p) / (r_a + r_p);
            alt_km = peri_alt_km;
            inc_deg = cfg.eccentric_inc_range(1) + diff(cfg.eccentric_inc_range) * rand();
            omega_deg = cfg.eccentric_omega_range(1) + diff(cfg.eccentric_omega_range) * rand();
            Omega_deg = cfg.eccentric_RAAN_range(1) + diff(cfg.eccentric_RAAN_range) * rand();
        case 'polar'
            alt_km = alt_range(1) + diff(alt_range) * rand();
            is_circular = rand() < cfg.polar_circular_probability;
            if is_circular
                ecc = 0;
            else
                r_p = (params.R_moon_km + alt_km) * 1000;
                apo_offset = cfg.polar_apo_offset_km(1) + diff(cfg.polar_apo_offset_km) * rand();
                r_a = r_p + apo_offset * 1000;
                ecc = (r_a - r_p) / (r_a + r_p);
            end
            inc_deg = cfg.polar_inc_range(1) + diff(cfg.polar_inc_range) * rand();
            omega_deg = cfg.polar_omega_options(randi(numel(cfg.polar_omega_options)));
            Omega_deg = cfg.polar_RAAN_range(1) + diff(cfg.polar_RAAN_range) * rand();
        case 'general'
            alt_km = cfg.general_alt_range(1) + diff(cfg.general_alt_range) * rand();
            ecc = cfg.general_ecc_range(1) + diff(cfg.general_ecc_range) * rand();
            inc_deg = cfg.general_inc_range(1) + diff(cfg.general_inc_range) * rand();
            omega_deg = cfg.general_omega_range(1) + diff(cfg.general_omega_range) * rand();
            Omega_deg = cfg.general_RAAN_range(1) + diff(cfg.general_RAAN_range) * rand();
        case 'higheccentric'
            peri_alt_km = cfg.high_ecc_peri_range(1) + diff(cfg.high_ecc_peri_range) * rand();
            apo_alt_km = cfg.high_ecc_apo_range(1) + diff(cfg.high_ecc_apo_range) * rand();
            if apo_alt_km <= peri_alt_km, apo_alt_km = peri_alt_km + 1; end
            r_p = (params.R_moon_km + peri_alt_km) * 1000;
            r_a = (params.R_moon_km + apo_alt_km) * 1000;
            ecc = (r_a - r_p) / (r_a + r_p);
            alt_km = peri_alt_km;
            inc_deg = cfg.high_ecc_inc_range(1) + diff(cfg.high_ecc_inc_range) * rand();
            omega_deg = cfg.high_ecc_omega_range(1) + diff(cfg.high_ecc_omega_range) * rand();
            Omega_deg = cfg.high_ecc_RAAN_range(1) + diff(cfg.high_ecc_RAAN_range) * rand();
        case 'frozen'
            % Literature-based frozen orbit generation
            % Based on Folta & Quinn (2006), Elipe & Lara (2003), Nie & Gurfil (2018)

            % Select frozen orbit family based on weights
            family_names = fieldnames(cfg.frozen_families);
            cum_weights = cumsum(cfg.frozen_family_weights);
            
            % Debug: Show family selection details (only for first orbit to avoid spam)
            if i <= 3 && strcmp(family_type, 'Frozen')  % Show first 3 selections for debugging
                %fprintf('Debug family selection (orbit %d):\n', i);
                %fprintf('  Family names: '); fprintf('%s ', family_names{:}); fprintf('\n');
                %fprintf('  Weights: '); fprintf('%.3f ', cfg.frozen_family_weights); fprintf('\n');
                %fprintf('  Cum weights: '); fprintf('%.3f ', cum_weights); fprintf('\n');
            end
            
            r = rand();
            selected_idx = find(r <= cum_weights, 1, 'first');
            selected_family = family_names{selected_idx};
            
            if i <= 3 && strcmp(family_type, 'Frozen')  % Show first 3 selections for debugging
                %fprintf('  Random r=%.4f -> Selected index %d -> Family %s\n\n', r, selected_idx, selected_family);
            end
            
            family_params = cfg.frozen_families.(selected_family);

            % Set altitude - use family-specific range if available
            if isfield(family_params, 'alt_range_km')
                alt_km = family_params.alt_range_km(1) + diff(family_params.alt_range_km) * rand();
            else
                alt_km = alt_range(1) + diff(alt_range) * rand();
            end

            % Set inclination (fixed for each family, with tolerance-based variation)
            if isfield(family_params, 'inc_tolerance_deg')
                inc_tolerance = family_params.inc_tolerance_deg;
            else
                inc_tolerance = 2.0; % Default tolerance
            end
            inc_deg = family_params.inc_deg + inc_tolerance * (rand() - 0.5);

            % Set eccentricity based on family
            % All families now use ecc_range for consistency and stability
            if isfield(family_params, 'ecc_range')
                ecc = family_params.ecc_range(1) + diff(family_params.ecc_range) * rand();
            else
                % Handle families without ecc_range (e.g., polar)
                if strcmp(selected_family, 'polar')
                    % For polar frozen, use small eccentricity range
                    ecc = 0.001 + 0.019 * rand(); % [0.001, 0.02] similar to magic
                else
                    ecc = 0.01; % default fallback
                end
            end

            % Argument of perigee: should librate around ±90°
            % Use family-specific libration amplitude if available
            if isfield(family_params, 'omega_libration_deg')
                libration_amp = family_params.omega_libration_deg;
            else
                libration_amp = 5.0; % Default libration amplitude
            end
            
            if rand() < 0.5
                omega_deg = 90 + libration_amp * (rand() - 0.5);
            else
                omega_deg = 270 + libration_amp * (rand() - 0.5);
            end

            % RAAN: full range
            Omega_deg = cfg.frozen_RAAN_range(1) + diff(cfg.frozen_RAAN_range) * rand();

            % Store family information for tracking
            family_name_for_storage = sprintf('Frozen_%s', selected_family);
    end
    
    nu_deg = 0; % Start at periapsis

    % --- Conversion to Cartesian State in MCI Frame (SI Units) ---
    r_p_m = (params.R_moon_km + alt_km) * 1000;
    
    % Safeguard for high eccentricity orbits
    if ecc >= 0.99
        % For near-parabolic orbits, use a different approach
        warning('High eccentricity detected (e=%.4f), adjusting semi-major axis calculation', ecc);
        a_m = r_p_m / (1 - 0.95); % Use a safer eccentricity for calculation
    else
        a_m = r_p_m / (1 - ecc);
    end

    %fprintf('Orbit generation: alt=%.1f km, ecc=%.4f, a=%.1f km\n', alt_km, ecc, a_m/1000);

    [r_vec, v_vec] = keplerian_to_cartesian_MCI(a_m, ecc, ...
        deg2rad(inc_deg), deg2rad(omega_deg), deg2rad(Omega_deg), deg2rad(nu_deg), params.GM_moon);

    % Store results
    if strcmpi(family_type, 'Frozen') && exist('family_name_for_storage', 'var')
        Families{i} = family_name_for_storage;
    else
        Families{i} = family_type;
    end
    Altitudes(i) = alt_km;
    Inclinations(i) = inc_deg;
    Eccentricities(i) = ecc;
    ArgPeriapsis(i) = omega_deg;
    RAANs(i) = Omega_deg;
    TrueAnomalies(i) = nu_deg;
    X(i) = r_vec(1); Y(i) = r_vec(2); Z(i) = r_vec(3);
    VX(i) = v_vec(1); VY(i) = v_vec(2); VZ(i) = v_vec(3);
    
    % Validate the generated orbit
    E_check = (norm(v_vec)^2 / 2) - params.GM_moon / norm(r_vec);
    if E_check >= 0
        warning('Generated orbit is not bound! Energy = %.6e J/kg', E_check);
    end
    if norm(r_vec) < params.R_moon
        warning('Generated orbit intersects Moon surface! r = %.1f km', norm(r_vec)/1000);
    end
end

primaries_table = table(IDs, Families, Altitudes, Inclinations, Eccentricities, ...
    ArgPeriapsis, RAANs, TrueAnomalies, X, Y, Z, VX, VY, VZ, ...
    'VariableNames', {'ID', 'Family', 'Altitude_km', 'Inclination_deg', 'Eccentricity', ...
    'ArgPeriapsis_deg', 'RAAN_deg', 'TrueAnomaly_deg', 'x', 'y', 'z', 'vx', 'vy', 'vz'});
end

function [r_vec, v_vec] = keplerian_to_cartesian_MCI(a, e, inc, omega, Omega, nu, mu_central)
% Convert Keplerian elements to Cartesian state in an Inertial Frame (SI units)
p = a * (1 - e^2);
if p < 1e-6, p = a; end % Handle circular case

r_mag = p / (1 + e * cos(nu));
r_peri = [r_mag * cos(nu); r_mag * sin(nu); 0];
v_peri = [-sqrt(mu_central / p) * sin(nu); sqrt(mu_central / p) * (e + cos(nu)); 0];

R_Omega = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1];
R_i = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
R_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];

R_peri_to_I = (R_Omega * R_i * R_omega)';
r_vec = R_peri_to_I * r_peri;
v_vec = R_peri_to_I * v_peri;
end

function cfg_out = apply_overrides(cfg_in, ov)
cfg_out = cfg_in;
if isempty(ov), return; end
f = fieldnames(ov);
for k=1:numel(f)
    cfg_out.(f{k}) = ov.(f{k});
end
end
