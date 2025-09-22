% Canonical copy from LLO_population/dynamicsSTM_perturbed_2bp.m
% DO NOT MODIFY HERE.
function d_aug_state = dynamicsSTM_perturbed_2bp(t, aug_state, params)
%DYNAMICSSTM_PERTURBED_2BP Dynamics for state and STM propagation.
%
% This function computes the derivative of an augmented state vector, which
% includes the 6-element spacecraft state and the 36-element state
% transition matrix (STM). It includes perturbations from lunar harmonics
% (J2, J3, C22) and third-body gravity from the Earth.

% Input validation
if ~isscalar(t)
    error('t must be a scalar');
end
if ~isvector(aug_state) || length(aug_state) ~= 42
    error('aug_state must be a 42-element vector');
end
if ~isstruct(params)
    error('params must be a struct');
end

% --- 1. Unpack State and STM ---
state = aug_state(1:6);
STM = reshape(aug_state(7:42), 6, 6);

% --- 2. Compute Acceleration from Verified RHS Function ---
% The state derivative (velocity and acceleration) is calculated by the 
% trusted two_body_perturbed_rhs function.
d_state = two_body_perturbed_rhs(t, state, params);

% --- 3. Build the Jacobian Matrix A(t) ---
% The Jacobian A(t) for the variational equations d(STM)/dt = A(t)*STM
% has the structure:
% A = [ 0   I ]
%     [ G   0 ]
% where G is the 3x3 gravity gradient matrix d(acceleration)/d(position).

A = zeros(6, 6);
A(1:3, 4:6) = eye(3);

% Compute gravity gradient numerically from the verified total acceleration model
% This ensures consistency with J2 and all other perturbations used in the RHS.
r_sc_I = state(1:3);
v_sc_I = state(4:6);
r_norm = norm(r_sc_I);
% Finite-difference step scaled to radius to balance truncation and rounding errors
h = max(1e-3, 1e-6 * r_norm); % meters
G = zeros(3,3);
e = eye(3);
for i = 1:3
    state_p = [r_sc_I + h * e(:,i); v_sc_I];
    state_m = [r_sc_I - h * e(:,i); v_sc_I];
    d_state_p = two_body_perturbed_rhs(t, state_p, params);
    d_state_m = two_body_perturbed_rhs(t, state_m, params);
    a_p = d_state_p(4:6);
    a_m = d_state_m(4:6);
    G(:, i) = (a_p - a_m) / (2*h);
end
A(4:6, 1:3) = G;

% --- 4. Propagate the STM ---
% Reference: Vallado, Eq. 8-72
STM_dot = A * STM;

% --- 5. Assemble the Augmented State Derivative ---
d_aug_state = [d_state; reshape(STM_dot, 36, 1)];
end

