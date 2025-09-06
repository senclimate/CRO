function par_cell = par_load(data_name, ro_name)
% PAR_LOAD  Load CRO parameters by (data_name, ro_name) from a fixed .mat file.
% - If data_name ends with "-all", it loads all "<base>-<number>" entries for the given ro_name
%   and concatenates them column-wise into a 16xN cell (sorted by <number> asc).
% - Otherwise, it requires exactly one match and returns the original 16x1 cell.
%
% Usage:
%   p  = par_load("ORAS5","Linear-White-Additive");                 % returns 16x1 cell
%   pa = par_load("CMIP6-historical-all","Linear-White-Additive");  % returns 16xN cell (e.g., N=48)
%
% Notes:
%   Expects CRO_parlib_v0.0.mat with variables: par, tag_data_type, tag_RO_type.
%   Each S.par{i,j} is a 16x1 cell whose entries can be scalars/vectors/matrices.

    FNAME = "CRO_parlib_v0.0.mat";  % Change this path/name here if needed.

    % Load required variables from the .mat file
    S = load(FNAME, 'par', 'tag_data_type', 'tag_RO_type');

    data_name = string(data_name);
    ro_name   = string(ro_name);

    % ---- Case 1: "-all" suffix -> load all numbered variants ----
    if endsWith(data_name, "-all")
        base = erase(data_name, "-all");  % prefix before the numbering

        names = string(S.tag_data_type);
        ros   = string(S.tag_RO_type);

        % Build case-insensitive regex: ^<base>-(\d+)$
        pat = "^" + regexptranslate('escape', base) + "-(\d+)$";

        % Find names matching "<base>-<number>" AND exact ro_name
        m = ~cellfun('isempty', regexp(names, pat, 'once', 'ignorecase')) & (ros == ro_name);
        idxs = find(m);

        if isempty(idxs)
            error('par_load:NotFoundAll', ...
                  'No numbered records found for base="%s", RO_type="%s".', base, ro_name);
        end

        % Extract numeric suffixes and sort ascending
        nums = nan(numel(idxs),1);
        for k = 1:numel(idxs)
            tok = regexp(names(idxs(k)), pat, 'tokens', 'once', 'ignorecase');
            nums(k) = str2double(tok{1});
        end

        [~, order] = sort(nums, 'ascend');
        idxs = idxs(order);

        % Concatenate each 16x1 cell column-wise into 16xN
        par_cell = [S.par{idxs}];  % brace expansion concatenates contents

        return
    end

    % ---- Case 2: exact single match required (no "-all") ----
    % Use exact (case-sensitive) match as you requested "exactly one"
    idx = strcmp(string(S.tag_data_type), data_name) & ...
          strcmp(string(S.tag_RO_type),   ro_name);

    match_count = nnz(idx);
    if match_count == 0
        error('par_load:NotFound', ...
              'No record found for data="%s", RO_type="%s".', data_name, ro_name);
    elseif match_count > 1
        error('par_load:MultipleMatches', ...
              'Multiple records found for data="%s", RO_type="%s".', data_name, ro_name);
    end

    [i, j]  = find(idx, 1, 'first');
    par_cell = S.par{i, j};  % returns the original 16x1 cell
end
