function prim = load_database_primaries(dbRoot)
%LOAD_DATABASE_PRIMARIES Load periodic-orbit primaries from Database_Orbits.
% Usage:
%   prim = load_database_primaries('C:\Users\<you>\Desktop\cislunar-test-case-creator-main\Database_Orbits');
%
% Scans known orbit-family subfolders, reads CSVs, and aggregates:
%   - CR3BP initial conditions x0 = [x y z vx vy vz] (normalized units)
%   - Jacobi constant C (LU^2/TU^2)
%   - Period T (TU) and Period_days (if present)
%   - Mass ratio mu (if present)
% Families aggregated (fields in output struct):
%   DRO, L1_halo, L2_halo, L1_lyapunov, L2_lyapunov
%
% Output prim.<family> has fields:
%   .name       family name
%   .states     6xN matrix of initial states (CR3BP, normalized)
%   .C          1xN Jacobi constants
%   .T          1xN periods in TU (NaN if unavailable)
%   .T_days     1xN periods in days (NaN if unavailable)
%   .mu         1xN mass ratio values (NaN if unavailable)
%   .ids        1xN cell array of IDs (if available)
%   .sources    1xN cell array of source filenames for each row
%   .count      number of rows aggregated
%   .notes      info text
%
% Notes:
% - Header parsing is robust to variants (e.g., 'x0 (LU)', 'vx0 (LU/TU)')
% - Only CSV files are considered; MAT support can be added if needed.

if nargin < 1 || ~(ischar(dbRoot) || isstring(dbRoot))
    error('load_database_primaries:arg', 'Provide dbRoot folder path as char or string.');
end
dbRoot = char(dbRoot);
if exist(dbRoot, 'dir') ~= 7
    error('load_database_primaries:nodir', 'Folder not found: %s', dbRoot);
end

families = {
    'DRO',        {'DROs','DRO'};
    'L1_halo',    {'L1_Halo','L1_halo'};
    'L2_halo',    {'L2_Halo','L2_halo'};
    'L1_lyapunov',{'L1_Lyapunov','L1_lyapunov'};
    'L2_lyapunov',{'L2_Lyapunov','L2_lyapunov'}
};

prim = struct();
prim.root = dbRoot;
prim.families = families(:,1)';

for i = 1:size(families,1)
    famName = families{i,1};
    subCandidates = families{i,2};
    folder = '';
    for k = 1:numel(subCandidates)
        p = fullfile(dbRoot, subCandidates{k});
        if exist(p, 'dir') == 7
            folder = p; break
        end
    end
    fam = struct('name', famName, 'states', [], 'C', [], 'T', [], 'T_days', [], ...
                 'mu', [], 'ids', {{}}, 'sources', {{}}, 'count', 0, 'notes', "");
    if isempty(folder)
        fam.notes = "No matching subfolder found";
        prim.(famName) = fam; %#ok<*STRNU>
        continue
    end

    csvFiles = dir(fullfile(folder, '**', '*.csv'));
    if isempty(csvFiles)
        fam.notes = "No CSV files found";
        prim.(famName) = fam;
        continue
    end

    states = [];
    C = []; T = []; T_days = []; mu = [];
    ids = {}; sources = {};

    for c = 1:numel(csvFiles)
        fp = fullfile(csvFiles(c).folder, csvFiles(c).name);
        try
            Ttab = readtable(fp, 'VariableNamingRule','preserve');
        catch ME
            warning('load_database_primaries:read', 'Failed reading %s (%s)', fp, ME.identifier);
            continue
        end
        if height(Ttab) == 0, continue; end

        idx = resolve_columns_extended(Ttab);
        % Extract rows
        Xi = [Ttab{:, idx.x},  Ttab{:, idx.y},  Ttab{:, idx.z}, ...
              Ttab{:, idx.vx}, Ttab{:, idx.vy}, Ttab{:, idx.vz}];
        states = [states; Xi]; %#ok<AGROW>

        C  = [C;  fetchOrNaN(Ttab, idx.C)]; %#ok<AGROW>
        Ti = fetchOrNaN(Ttab, idx.T);
        Td = fetchOrNaN(Ttab, idx.T_days);
        mui= fetchOrNaN(Ttab, idx.mu);
        T   = [T; Ti]; %#ok<AGROW>
        T_days = [T_days; Td]; %#ok<AGROW>
        mu  = [mu; mui]; %#ok<AGROW>

        if ~isnan(idx.id)
            idVals = Ttab{:, idx.id};
            if isnumeric(idVals)
                idVals = arrayfun(@(v) sprintf('%g', v), idVals, 'UniformOutput', false);
            elseif isstring(idVals)
                idVals = cellstr(idVals);
            end
        else
            idVals = repmat({''}, height(Ttab), 1);
        end
        ids = [ids; idVals]; %#ok<AGROW>
        sources = [sources; repmat({fp}, height(Ttab), 1)]; %#ok<AGROW>
    end

    % Transpose into 6xN and row vectors
    if ~isempty(states)
        fam.states  = states';           % 6xN
        fam.C       = C.';               % 1xN
        fam.T       = T.';               % 1xN
        fam.T_days  = T_days.';          % 1xN
        fam.mu      = mu.';              % 1xN
        fam.ids     = ids.';             % 1xN cell
        fam.sources = sources.';         % 1xN cell
        fam.count   = size(states,1);
        fam.notes   = sprintf('Loaded %d rows from %d CSV file(s)', fam.count, numel(csvFiles));
    else
        fam.notes   = sprintf('No rows parsed from %d CSV file(s)', numel(csvFiles));
    end

    prim.(famName) = fam;
end

end

function idx = resolve_columns_extended(T)
% Build robust mapping from table headers to indices for required/optional columns.
    vnOrig = string(T.Properties.VariableNames);
    norm = @(s) regexprep(lower(s), '[^a-z0-9]', '');
    vn = norm(vnOrig);
    findFirst = @(patterns) local_find(vn, patterns, norm);

    idx = struct();
    % Position & velocity (required)
    idx.x  = findFirst({"x0","x","posx","xlu"});
    idx.y  = findFirst({"y0","y","posy","ylu"});
    idx.z  = findFirst({"z0","z","posz","zlu"});
    idx.vx = findFirst({"vx0","vx","velx","v_x","vxlu","lutuvx"});
    idx.vy = findFirst({"vy0","vy","vely","v_y","vylu","lutuvy"});
    idx.vz = findFirst({"vz0","vz","velz","v_z","vzlu","lutuvz"});

    req = [idx.x, idx.y, idx.z, idx.vx, idx.vy, idx.vz];
    if any(isnan(req))
        err = "";
        names = {"x","y","z","vx","vy","vz"};
        for k=1:numel(req)
            if isnan(req(k)), err = err + names{k} + " "; end
        end
        error('load_database_primaries:missingColumns', 'Missing required columns: %s', err);
    end

    % Optional columns
    idx.C      = findFirst({"jacobi","j","c","jacobiconstant"});
    idx.T      = findFirst({"periodtu","period","t","tperiod"});
    idx.T_days = findFirst({"perioddays","days","tday","periodindays"});
    idx.mu     = findFirst({"massratio","mu"});
    idx.id     = findFirst({"id","index","name"});
    idx.stab   = findFirst({"stabilityindex","stability","sindex"}); %#ok<NASGU>
end

function j = local_find(vn, patterns, normfun)
% Helper: find first header matching any pattern (prefix match ok).
    for p = 1:numel(patterns)
        cand = normfun(string(patterns{p}));
        j = find(startsWith(vn, cand), 1);
        if ~isempty(j), return; end
    end
    j = NaN;
end

function v = fetchOrNaN(Ttab, idx)
    if isnan(idx), v = nan(height(Ttab),1); else, v = Ttab{:, idx}; end
end
