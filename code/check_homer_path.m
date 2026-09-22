function check_homer_path()
% CHECK_HOMER_PATH  Ensure Homer3 is on the MATLAB path.
%
% The submitted version hard-coded a local Homer3 folder. It now resolves the
% location through setup_paths_config.m so the repository is portable.
% Homer3 v1.87.0 was used for the published analyses.

if ~isempty(which('SnirfLoad'))
    return
end

P = setup_paths();
if isfield(P, 'homer3') && ~isempty(P.homer3) && isfolder(P.homer3)
    addpath(genpath(P.homer3));
end

assert(~isempty(which('SnirfLoad')), ...
    ['Homer3 not found. Copy setup_paths_config_example.m to ' ...
     'setup_paths_config.m and set the homer3 field. Homer3 v1.87.0 was used ' ...
     'for the published analyses: https://github.com/BUNPC/Homer3']);
end

% EOF
