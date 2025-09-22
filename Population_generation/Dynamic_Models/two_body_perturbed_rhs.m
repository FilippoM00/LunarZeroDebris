function d_state = two_body_perturbed_rhs(t, state, params)
% TWO_BODY_PERTURBED_RHS  Moon-centered inertial dynamics with J2 and Earth 3rd-body.
% State: [x;y;z;vx;vy;vz] in meters and m/s. Parameters in SI units.
%
% Aligns with PoC_Analysis LLO dynamics for consistency:
% - Central Moon gravity
% - Earth third-body differential acceleration (time-varying Earth position)
% - Lunar J2 applied in Moon-fixed frame and rotated back to inertial

% Unpack state
r_I = state(1:3);
v_I = state(4:6);

%% Central Moon gravity
mu_m = params.GM_moon;
rnorm = norm(r_I);
a_moon_I = -mu_m * r_I / rnorm^3;

%% Earth third-body (time-varying on circular orbit in MCI XY-plane)
mu_e = params.GM_earth;
D_em = params.D_em;
n_e  = params.n_earth; % mean motion of Earth around Moon

theta_e = n_e * t + pi; % Earth starts at -X at t=0 (matches PoC)
r_me_I  = [D_em * cos(theta_e); D_em * sin(theta_e); 0];  % Moon->Earth
r_se_I  = r_me_I - r_I;                                   % SC->Earth

% Differential acceleration (Curtis eqn.):
% a3b = mu_e * ( r_se/|r_se|^3 - r_me/|r_me|^3 )
a_3b_I = mu_e * ( r_se_I / norm(r_se_I)^3 - r_me_I / norm(r_me_I)^3 );

%% Lunar J2 (axisymmetric; computed in Moon-fixed frame, then rotated to inertial)
J2   = getfield_safe(params, 'J2', 0);
R    = getfield_safe(params, 'R_moon', 1737.4e3);
wm   = getfield_safe(params, 'w_moon', 2.6617e-6);

if J2 ~= 0
	% Inertial->Fixed rotation (Moon rotates about +Z)
	theta_m = wm * t;
	R_I2F = [ cos(theta_m)  sin(theta_m)  0; 
			 -sin(theta_m)  cos(theta_m)  0;
			  0             0             1];
	r_F = R_I2F * r_I;
	x = r_F(1); y = r_F(2); z = r_F(3);
	r2 = x*x + y*y + z*z;
	r5 = r2^(2.5);
	z2 = z*z;
	f = -1.5 * J2 * mu_m * R^2 / r5; % factor with sign chosen to match standard form
	ax_F = f * x * (1 - 5*z2/r2);
	ay_F = f * y * (1 - 5*z2/r2);
	az_F = f * z * (3 - 5*z2/r2);
	a_J2_I = R_I2F' * [ax_F; ay_F; az_F];
else
	a_J2_I = [0;0;0];
end

%% Total acceleration and derivative
a_total_I = a_moon_I + a_3b_I + a_J2_I;
d_state = [v_I; a_total_I];
end

function val = getfield_safe(s, name, default)
if isfield(s, name)
	val = s.(name);
else
	val = default;
end
end
