function P = setup_paths_config_example()
% Locations of the external toolboxes, on THIS machine.
% Copy this file to setup_paths_config.m and edit the paths.
% setup_paths_config.m is gitignored, so local paths are never committed.
%
% Leave a field empty ('') if you do not have that toolbox; only Homer3, the
% Brain Connectivity Toolbox and the Statistics and Machine Learning Toolbox
% are needed for the revision-01 pipeline.

P = struct( ...
    'homer3',        '', ...   % e.g. 'C:\toolboxes\Homer3'
    'bct',           '', ...   % e.g. 'C:\toolboxes\BCT'
    'shadedErrorBar','', ...   % e.g. 'C:\toolboxes\shadedErrorBar'
    'qtnirs',        '', ...   % e.g. 'C:\toolboxes\qt-nirs'
    'nbs',           '', ...   % e.g. 'C:\toolboxes\NBS1.2'
    'fieldtrip',     '', ...   % e.g. 'C:\toolboxes\fieldtrip'
    'circularGraph', '');      % optional, one plotting script only
end

% EOF
