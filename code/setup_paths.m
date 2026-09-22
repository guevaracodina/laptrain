function P = setup_paths(cfgFile)
% SETUP_PATHS  Add the external toolboxes this project depends on.
%
%   setup_paths            adds every configured toolbox to the MATLAB path
%   P = setup_paths        also returns the resolved locations
%   P = setup_paths(file)  reads locations from a different config file
%
% HOW TO USE
%   1. Copy setup_paths_config_example.m to setup_paths_config.m
%   2. Edit it so each field points at your own copy of the toolbox
%   3. Call setup_paths once per MATLAB session, before running anything else
%
% setup_paths_config.m is listed in .gitignore, so your local paths are never
% committed. If it is absent, this function falls back to whatever is already
% on the MATLAB path and reports what is missing rather than failing silently.
%
% REQUIRED, for the revision-01 pipeline
%   homer3          Homer3 v1.87            https://github.com/BUNPC/Homer3
%   bct             Brain Connectivity Tbx  https://sites.google.com/site/bctnet/
% REQUIRED for specific scripts only
%   shadedErrorBar  fc_network_measures_*   https://github.com/raacampbell/shadedErrorBar
%   qtnirs          qt_nirs_script          https://github.com/lpollonini/qt-nirs
%   nbs             NBS_script (GUI)        https://www.nitrc.org/projects/nbs/
%   fieldtrip       code/BIDS/*             https://www.fieldtriptoolbox.org/
%   circularGraph   fc_network_measures_plot_script only
%
% The Statistics and Machine Learning Toolbox is also required (fitlme, fitrm,
% multcompare, ttest, ranksum). It ships with MATLAB and needs no path entry.

if nargin < 1 || isempty(cfgFile), cfgFile = 'setup_paths_config'; end

if exist(cfgFile, 'file') == 2
    P = feval(cfgFile);
else
    warning('setup_paths:noConfig', ...
        ['%s.m not found. Copy setup_paths_config_example.m to %s.m and edit it.\n' ...
         'Falling back to toolboxes already on the MATLAB path.'], cfgFile, cfgFile);
    P = struct();
end

% Add each configured folder, recursively where the toolbox needs it
recursive = {'homer3', 'bct', 'qtnirs', 'nbs', 'fieldtrip'};
fn = fieldnames(P);
for k = 1:numel(fn)
    d = P.(fn{k});
    if isempty(d), continue; end
    if ~isfolder(d)
        warning('setup_paths:missingFolder', '%s: not a folder: %s', fn{k}, d);
        continue
    end
    if ismember(fn{k}, recursive), addpath(genpath(d)); else, addpath(d); end
end

% Report what is actually reachable, so a missing toolbox surfaces now and not
% forty minutes into an analysis
checks = { 'Homer3',                      'SnirfLoad'
           'Brain Connectivity Toolbox',  'threshold_proportional'
           'shadedErrorBar',              'shadedErrorBar'
           'QT-NIRS',                     'qtnirs'
           'Statistics and ML Toolbox',   'fitlme' };
fprintf('\nToolbox availability\n');
missingRequired = {};
for k = 1:size(checks, 1)
    ok = ~isempty(which(checks{k,2}));
    if ok, status = 'found'; else, status = 'NOT FOUND'; end
    fprintf('  %-28s %s\n', checks{k,1}, status);
    if ~ok && ismember(checks{k,1}, {'Homer3','Brain Connectivity Toolbox', ...
                                     'Statistics and ML Toolbox'})
        missingRequired{end+1} = checks{k,1}; %#ok<AGROW>
    end
end
if ~isempty(missingRequired)
    warning('setup_paths:missingRequired', ...
            'Required toolbox(es) not on the path: %s', strjoin(missingRequired, ', '));
end
fprintf('\n');
end

% EOF
