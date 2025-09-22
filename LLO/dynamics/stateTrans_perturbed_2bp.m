% Canonical copy from LLO_population/stateTrans_perturbed_2bp.m
% NOTE: Modified to use standard ode113 output capture (no OutputFcn).

function [state_f, STM_f] = stateTrans_perturbed_2bp(state_i, tf, p, options, varargin)
% stateTrans_perturbed_2bp Propagates a state vector and the associated
% state transition matrix (STM) using a high-fidelity perturbed two-body model.
%
%   [state_f, STM_f] = stateTrans_perturbed_2bp(state_i, tf, p, options)
%
%   Inputs:
%   - state_i: 6x1 initial state vector [rx; ry; rz; vx; vy; vz] in the
%              Moon-Centered Inertial (MCI) frame (m, m/s).
%   - tf: Final propagation time (s). Can be positive (forward) or negative (backward).
%   - p: Struct of physical and simulation parameters, loaded from LLO_defaults.
%   - options: odeset options for the numerical integrator (e.g., from odeset).
%
%   Outputs:
%   - state_f: 6x1 final state vector at time tf.
%   - STM_f: 6x6 final state transition matrix at time tf.
%
%   This function initializes the 6x6 STM as the identity matrix, reshapes it
%   into a 36x1 vector, and appends it to the state vector to create a
%   42-element augmented state. It then calls ode113 with the
%   dynamicsSTM_perturbed_2bp function to integrate the full system.

% Input validation
if ~isvector(state_i) || length(state_i) ~= 6
    error('state_i must be a 6-element vector');
end
if ~isscalar(tf)
    error('tf must be a scalar');
end
if ~isstruct(p)
    error('p must be a struct');
end
if ~isstruct(options)
    error('options must be a struct from odeset');
end

% Optional absolute time offset t0 so dynamics evaluate at t + t0
t0 = 0;
if ~isempty(varargin)
	t0 = varargin{1};
end

% Handle zero-time propagation quickly
if tf == 0
    state_f = state_i(:);
    STM_f = eye(6);
    return;
end

% Define the time span for integration (relative). Supports forward/backward.
tspan = [0, tf];

% Initialize the 6x6 State Transition Matrix (STM) as the identity matrix
STM_i = eye(6);

% Reshape the STM into a 36x1 column vector
STM_i_vec = reshape(STM_i, 36, 1);

% Create the 42x1 augmented initial state vector [state; stm]
augmented_state_i = [state_i; STM_i_vec];

% Propagate the augmented state using ode113 with absolute-time offset
% Capture the full output and extract only the final state
[~, Y] = ode113(@(t,y) dynamicsSTM_perturbed_2bp(t + t0, y, p), tspan, augmented_state_i, options);
final_state = Y(end, :)';
state_f = final_state(1:6);
STM_f = reshape(final_state(7:42), 6, 6);

end


