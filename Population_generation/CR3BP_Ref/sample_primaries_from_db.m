function S = sample_primaries_from_db(dbRoot, opts)
%SAMPLE_PRIMARIES_FROM_DB Build a sampled primaries set from Database_Orbits.
% Usage:
%   S = sample_primaries_from_db('C:\...\Database_Orbits');
%   S = sample_primaries_from_db(dbRoot, struct('seed',42));
%
% Returns a struct with fields:
%   .DRO, .L1_halo, .L2_halo, .L1_lyapunov, .L2_lyapunov, each containing:
%       states (6xN), C, T, T_days, mu, ids, sources, count
%   .counts_requested  map of requested counts
%   .counts_actual     map of actual counts
%   .notes             info
%
% Policy:
% - DRO: 1000 samples (without replacement if enough; with replacement otherwise)
% - Halos L1/L2: 1000 each, balanced across Northern/Southern if possible
% - Lyapunov L1/L2: 1000 each
% - Entries without a valid TU period (T) are discarded during sampling and backfilled

if nargin < 2, opts = struct(); end
if ~isfield(opts,'seed'), opts.seed = 12345; end
if ~isfield(opts,'counts')
    opts.counts = struct('DRO',100,'L1_halo',100,'L2_halo',100,'L1_lyapunov',100,'L2_lyapunov',100);
end
if ~isfield(opts,'enforce_dro_altitude_cap'), opts.enforce_dro_altitude_cap = true; end
if ~isfield(opts,'dro_altitude_cap_km'), opts.dro_altitude_cap_km = 100000; end
if ~isfield(opts,'strict_no_replacement_for_dro'), opts.strict_no_replacement_for_dro = true; end
if ~isfield(opts,'fill_shortfall_with_replacement'), opts.fill_shortfall_with_replacement = true; end
% L1 halo z-cap (normalized CR3BP units): enforce |z0| <= l1_halo_z_cap
if ~isfield(opts,'enforce_l1_halo_z_cap'), opts.enforce_l1_halo_z_cap = true; end
if ~isfield(opts,'l1_halo_z_cap'), opts.l1_halo_z_cap = 0.5; end
% L2 Lyapunov min-x filter (normalized CR3BP units): enforce x0 >= l2_lyapunov_min_x
if ~isfield(opts,'enforce_l2_lyapunov_min_x'), opts.enforce_l2_lyapunov_min_x = true; end
if ~isfield(opts,'l2_lyapunov_min_x'), opts.l2_lyapunov_min_x = 1.041; end
% L1 Lyapunov min-x filter (normalized CR3BP units): enforce x0 >= l1_lyapunov_min_x
if ~isfield(opts,'enforce_l1_lyapunov_min_x'), opts.enforce_l1_lyapunov_min_x = true; end
if ~isfield(opts,'l1_lyapunov_min_x'), opts.l1_lyapunov_min_x = 0.811; end

rng(opts.seed);
prim = load_database_primaries(dbRoot);

S = struct();
S.notes = sprintf('Sampled on %s with seed=%d', datestr(now), opts.seed);
S.counts_requested = opts.counts;

famList = {'DRO','L1_halo','L2_halo','L1_lyapunov','L2_lyapunov'};
for i=1:numel(famList)
    fam = famList{i};
    if ~isfield(prim, fam)
        warning('sample_primaries_from_db:missingFamily', 'Family not found in prim: %s', fam);
        S.(fam) = emptyFam(fam);
        continue
    end
    src = prim.(fam);
    % Require valid TU period
    valid = find(isfinite(src.T) & (src.T > 0));
    nReq = opts.counts.(fam);
    if strcmp(fam, 'L1_halo') || strcmp(fam, 'L2_halo')
        [idxN, idxS, ~] = split_halo_indices(src.sources);
        % Build the family-valid set (valid TU + optional z-cap for L1 halos)
        validFam = valid;
        if strcmp(fam, 'L1_halo') && opts.enforce_l1_halo_z_cap
            z0 = src.states(3, :);
            validZ = find(abs(z0) <= opts.l1_halo_z_cap);
            validFam = intersect(validFam, validZ, 'stable');
        end
        % Enforce family-valid set on N/S pools
        idxN = intersect(idxN, validFam, 'stable');
        idxS = intersect(idxS, validFam, 'stable');
        % Balanced pick N/S
        half = floor(nReq/2);
        pickN = pick_indices(idxN, numel(idxN), half);
        pickS = pick_indices(idxS, numel(idxS), nReq - half);
        idx = [pickN(:); pickS(:)];
        if numel(idx) < nReq
            % backfill from remaining valid pool if needed
            need = nReq - numel(idx);
            pool = setdiff(validFam, idx);
            idx = [idx; pick_indices(pool, numel(pool), need)];
            if strcmp(fam, 'L1_halo') && opts.enforce_l1_halo_z_cap && numel(idx) < nReq
                warning('sample_primaries_from_db:l1HaloShortfall', ...
                    'L1 halo shortfall: requested %d, available %d after filters (valid T + |z0|<=%.3f).', ...
                    nReq, numel(idx), opts.l1_halo_z_cap);
            end
        end
    elseif strcmp(fam, 'DRO') && opts.enforce_dro_altitude_cap
        % Filter DRO pool by max lunar altitude across 1 period
        pass = dro_within_altitude_cap(src, opts.dro_altitude_cap_km);
        pool = find(pass);
        nAvail = numel(pool);
        if opts.strict_no_replacement_for_dro
            if nAvail >= nReq
                idx = pick_without_replacement(pool, nReq);
            else
                % Underfilled: either fill remainder with replacement or leave short
                idx = pick_without_replacement(pool, nReq);
                if opts.fill_shortfall_with_replacement && nAvail > 0
                    need = nReq - numel(idx);
                    repl = pool(randi(nAvail, need, 1));
                    idx = [idx(:); repl(:)];
                elseif nAvail < nReq
                    warning('sample_primaries_from_db:droShortfall', ...
                        'DRO shortfall: requested %d, available %d after filters (valid T + altitude cap).', nReq, nAvail);
                end
            end
        else
            % Not strict: allow replacement to always reach nReq
            idx = pick_indices(pool, numel(pool), nReq);
        end
    elseif strcmp(fam, 'L2_lyapunov') && opts.enforce_l2_lyapunov_min_x
        % Filter L2 Lyapunov pool by minimum initial x-coordinate
        x0 = src.states(1, :);
        validFam = intersect(valid, find(x0 >= opts.l2_lyapunov_min_x), 'stable');
        idx = pick_indices(validFam, numel(validFam), nReq);
    elseif strcmp(fam, 'L1_lyapunov') && opts.enforce_l1_lyapunov_min_x
        % Filter L1 Lyapunov pool by minimum initial x-coordinate
        x0 = src.states(1, :);
        validFam = intersect(valid, find(x0 >= opts.l1_lyapunov_min_x), 'stable');
        idx = pick_indices(validFam, numel(validFam), nReq);
    else
        % Generic families: pick from valid TU set
        idx = pick_indices(valid, numel(valid), nReq);
    end
    S.(fam) = take_rows(src, idx);
end

% Summaries
for i=1:numel(famList)
    fam = famList{i}; S.counts_actual.(fam) = S.(fam).count; end

end

function pass = dro_within_altitude_cap(dro, cap_km)
% Return logical vector marking DRO rows with max lunar altitude <= cap_km
    pass = false(1, dro.count);
    if dro.count == 0, return; end
    % Constants for normalization
    try
        D_EM_km   = astroConstants(7);   % Earth-Moon distance [km]
        R_moon_km = astroConstants(30);  % Moon radius [km]
    catch
        % Fallback values if astroConstants not on path
        D_EM_km   = 384400;
        R_moon_km = 1738.0;
    end
    rmax_norm = (R_moon_km + cap_km) / D_EM_km;  % normalized LU
    % Ensure CR3BP ODE exists
    if exist('CR3BPunpert','file') ~= 2
        error('sample_primaries_from_db:missingCR3BP', 'CR3BPunpert.m not on MATLAB path');
    end
    % Integration options
    rel = 1e-8; abs = 1e-10; %#ok<NASGU>
    ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    for i=1:dro.count
        x0 = dro.states(:, i);
        mu = dro.mu(i);
        if isnan(mu) || mu <= 0
            mu = 0.012150668; % default Earth-Moon mass ratio
        end
        T = dro.T(i);
        % Reject entries without a valid TU period
        if isnan(T) || T <= 0
            pass(i) = false; continue
        end
    % Integrate over one period with tight tolerances (let solver pick steps)
    tspan = [0 T];
        try
            [~, X] = ode113(@(t,x) CR3BPunpert(t,x,mu), tspan, x0(:)', ode_opts);
        catch
            % Retry on integration error
            [~, X] = ode113(@(t,x) CR3BPunpert(t,x,mu), tspan, x0(:)', ode_opts);
        end
        % Dist to Moon in CR3BP (Moon at (1-mu,0,0))
        rx = X(:,1) - (1 - mu);
        ry = X(:,2);
        rz = X(:,3);
        rmoon = sqrt(rx.^2 + ry.^2 + rz.^2);
        if max(rmoon) <= rmax_norm
            pass(i) = true;
        end
    end
end

function idx = pick_without_replacement(pool, nPick)
    if isempty(pool), idx = []; return; end
    nAvail = numel(pool);
    k = min(nPick, nAvail);
    idx = pool(randperm(nAvail, k));
end

function F = emptyFam(name)
    F = struct('name',name,'states',zeros(6,0),'C',[],'T',[],'T_days',[], ...
               'mu',[],'ids',{{}},'sources',{{}},'count',0,'notes','empty');
end

function idx = pick_indices(pool, nAvail, nPick)
    if nAvail <= 0
        idx = []; return
    end
    if nAvail >= nPick
        perm = randperm(nAvail, nPick);
        idx = pool(perm);
    else
        % sample with replacement to reach nPick
        perm = randi(nAvail, nPick, 1);
        idx = pool(perm);
    end
end

function [idxN, idxS, idxRest] = split_halo_indices(sources)
% Classify halo rows into Northern/Southern using source filename hints.
    if isrow(sources), sources = sources'; end
    n = numel(sources);
    idxN = []; idxS = []; idxRest = [];
    for i=1:n
        s = lower(string(sources{i}));
        if contains(s, 'northern')
            idxN(end+1) = i; %#ok<AGROW>
        elseif contains(s, 'southern')
            idxS(end+1) = i; %#ok<AGROW>
        else
            idxRest(end+1) = i; %#ok<AGROW>
        end
    end
    % If we only have generic entries, push them into both pools evenly
    if isempty(idxN) && ~isempty(idxRest)
        idxN = idxRest; idxRest = [];
    end
    if isempty(idxS) && ~isempty(idxRest)
        idxS = idxRest; idxRest = [];
    end
end

function Y = take_rows(X, idx)
% Gather selected rows/cols from a family struct X
    if isempty(idx)
        Y = emptyFam(X.name); return
    end
    Y = struct();
    Y.name    = X.name;
    Y.states  = X.states(:, idx);
    Y.C       = X.C(idx);
    Y.T       = X.T(idx);
    Y.T_days  = X.T_days(idx);
    Y.mu      = X.mu(idx);
    Y.ids     = X.ids(idx);
    Y.sources = X.sources(idx);
    Y.count   = numel(idx);
    Y.notes   = sprintf('Sampled %d entries', Y.count);
end
