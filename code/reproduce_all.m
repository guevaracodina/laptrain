%% REPRODUCE_ALL  End-to-end reproduction of the revision-01 results.
%
% Neurophotonics NPH-260108-1. Run this from the 'code' folder.
%
% This is the ordered pipeline. Each stage checks that its inputs exist and
% stops with an explanatory message rather than failing deep inside an
% analysis, so a missing prerequisite is obvious immediately.
%
% STAGES
%   0  path setup and toolbox check
%   1  QT-NIRS signal-quality assessment      -> data/qtnirs/<cond>/QC_reportTable_<cond>.mat
%   2  condense those tables                  -> data/qc_summary.mat
%   3  motion-artifact summary                -> data/motion_artifacts_summary.mat
%   4  connectivity, metrics, statistics      -> data/rev01_results.{mat,txt}
%   5  channel pruning and run exclusions     -> printed, and returned in P
%   6  NBS inputs                             -> data/conn_Mat_HbO_NBS_<cond>.mat etc.
%
% Stages 1 and 3 are slow and only need running once; they read every
% derivative file. Stages 2 and 4-6 are comparatively quick. Set the flags
% below to skip stages already completed.
%
% After stage 6, run NBS itself in its own GUI on the saved arrays; NBS has no
% documented batch entry point and is therefore not scripted here.

clear; clc

doQTNIRS  = true;    % stage 1, slow
doMotion  = true;    % stage 3, slow
doRecompute = true;  % stage 4 step 0, re-derive connectivity from derivatives

dataDir = fullfile('..', 'data');

%% Stage 0 -----------------------------------------------------------------
fprintf('=== Stage 0: paths ===\n');
setup_paths;
addpath(fullfile(pwd, 'rev01'), '-begin');   % rev01 must shadow code/
w = which('conn_regress_global');
assert(contains(w, 'rev01'), ...
    ['rev01 is not first on the path. conn_regress_global resolves to:\n  %s\n' ...
     'Run addpath(fullfile(pwd,''rev01''), ''-begin'') and try again.'], w);
fprintf('Using corrected conn_regress_global from rev01.\n');

for c = {'lap01','lap02','lap03','resting01','resting02','resting03'}
    d = fullfile(dataDir, c{1}, 'derivatives', 'homer');
    assert(isfolder(d), ...
        ['Derivatives not found: %s\nExtract data/BIDSsource/%sderivatives.7z ' ...
         'into that folder. See data/README.md.'], d, c{1});
end
fprintf('Homer3 derivatives present for all six conditions.\n');

%% Stage 1 -----------------------------------------------------------------
if doQTNIRS
    fprintf('\n=== Stage 1: QT-NIRS quality assessment (slow) ===\n');
    assert(~isempty(which('qtnirs')), ...
        'QT-NIRS not on the path; set the qtnirs field in setup_paths_config.m');
    qt_nirs_script;
else
    fprintf('\n=== Stage 1 skipped ===\n');
end

%% Stage 2 -----------------------------------------------------------------
fprintf('\n=== Stage 2: condense the QT-NIRS report tables ===\n');
if isfile(fullfile(dataDir, 'qc_summary.mat'))
    fprintf('qc_summary.mat already present; skipping.\n');
else
    qc_summary_export;
end

%% Stage 3 -----------------------------------------------------------------
if doMotion && ~isfile(fullfile(dataDir, 'motion_artifacts_summary.mat'))
    fprintf('\n=== Stage 3: motion-artifact summary (slow) ===\n');
    load_motion_artifacts_script;
else
    fprintf('\n=== Stage 3 skipped or already present ===\n');
end

%% Stage 4 -----------------------------------------------------------------
fprintf('\n=== Stage 4: connectivity, network metrics and statistics ===\n');
fprintf('(rev01_run_all honours its own doRecompute flag; edit it there)\n');
rev01_run_all;

%% Stage 5 -----------------------------------------------------------------
fprintf('\n=== Stage 5: pruning and run exclusions ===\n');
P = rev01_pruning_report(dataDir, [], fullfile(dataDir, 'BIDSsource')); %#ok<NASGU>

%% Stage 6 -----------------------------------------------------------------
fprintf('\n=== Stage 6: NBS inputs ===\n');
NBS_script;

fprintf(['\nDone.\n' ...
         'Results: %s and %s\n' ...
         'Next: run NBS in its GUI on conn_Mat_HbO_NBS_<cond>.mat with the\n' ...
         'matching design_NBS_<cond>.mat and exchange_NBS_<cond>.mat.\n'], ...
        fullfile(dataDir,'rev01_results.mat'), fullfile(dataDir,'rev01_results.txt'));

% EOF
