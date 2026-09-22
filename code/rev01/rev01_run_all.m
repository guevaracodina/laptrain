%% REV01_RUN_ALL  Master script for Neurophotonics NPH-260108-1, revision 01.
%
% Runs every analysis required by the reviewers and writes two files:
%   ..\data\rev01_results.mat   compact numeric results (small, attachable)
%   ..\data\rev01_results.txt   the same results as a readable report
%
% HOW TO RUN
%   1. Extract the Homer3 derivatives so that the following paths exist:
%          ..\data\lap01\derivatives\homer\sub-01\nirs\...
%          ..\data\lap02\...  lap03\...  resting01\...  resting02\...  resting03\...
%      They are currently archived as
%          ..\data\BIDSsource\lap01derivatives.7z  and the five siblings.
%   2. From the 'code' folder:
%          addpath('..\paper\Neurophotonics\rev01', '-begin')
%          rev01_run_all
%
% The '-begin' matters: rev01 must precede 'code' on the path so that the
% corrected conn_regress_global.m is the one that is called.
%
% STEP 0 re-derives every connectivity matrix with the corrected
% short-separation channel indices [3 24]. It is the slow step, roughly the
% same cost as the original script1_conn_mat_group_level. Set doRecompute to
% false to skip it once it has been run.

clear; close all; format compact; clc

doRecompute = false;     % connectivity is unchanged by the binarization fix;
                         % set true only if the derivatives change
nPerms      = 5000;
seed        = 42;

dataDir = fullfile('..', 'data');
partDir = fullfile(dataDir, 'participants');

%% Path checks
if isempty(which('threshold_proportional'))
end
assert(~isempty(which('threshold_proportional')), ...
       'Brain Connectivity Toolbox not found on the path.');
assert(~isempty(which('get_channels_from_template')), ...
       'The project code folder is not on the path.');

% Confirm the corrected SSC regression is the one that will be called
w = which('conn_regress_global');
fprintf('Using conn_regress_global from:\n  %s\n', w);
src = fileread(w);
assert(contains(src, 'get_channels_from_template'), ...
      ['The OLD conn_regress_global.m is first on the path. ' ...
       'Run addpath(''..\\paper\\Neurophotonics\\rev01'', ''-begin'') first.']);

check_homer_path

%% Step 0. Re-derive connectivity with the corrected short-channel indices
if doRecompute
    fprintf('\n=== Step 0: recomputing connectivity matrices ===\n');
    conds = {'lap01', 'lap02', 'lap03', 'resting01', 'resting02', 'resting03'};
    tags  = {'beginLapConn', 'midLapConn', 'endLapConn', ...
             'beginRestingConn', 'midRestingConn', 'endRestingConn'};
    for iCond = 1:numel(conds)
        d = fullfile(dataDir, conds{iCond}, 'derivatives', 'homer');
        assert(isfolder(d), 'Derivatives not found: %s  (extract the .7z archives)', d);
        for iHb = 1:3
            conn_mat_group_level(d, iHb, tags{iCond});
        end
    end
    fprintf('Connectivity matrices recomputed.\n');
end

%% Step 1. Network metrics for both conditions and all three chromophores
fprintf('\n=== Step 1: network metrics ===\n');
conditions   = {'Lap', 'Resting'};
chromophores = {'HbO', 'HbR', 'HbT'};
NM = struct();
for iCond = 1:numel(conditions)
    for iHb = 1:numel(chromophores)
        key = [conditions{iCond} chromophores{iHb}];
        NM.(key) = rev01_network_measures(dataDir, conditions{iCond}, ...
                                          chromophores{iHb}, nPerms, seed);
    end
end

%% Step 2. Data-quality covariates
fprintf('\n=== Step 2: data-quality covariate analysis ===\n');
qcFile  = fullfile(dataDir, 'qc_summary.mat');
motFile = fullfile(dataDir, 'motion_artifacts_summary.mat');
assert(isfile(qcFile),  'Run qc_summary_export first; %s not found.', qcFile);
assert(isfile(motFile), '%s not found; run load_motion_artifacts_script.', motFile);

QC  = load(qcFile,  'SCIsubj', 'PSPsubj', 'goodChanPct', 'condNames');
MOT = load(motFile, 'pctMotion');

QA = struct();
for iCond = 1:numel(conditions)
    key = [conditions{iCond} 'HbO'];
    QA.(conditions{iCond}) = rev01_quality_covariate( ...
        NM.(key), QC.SCIsubj, QC.PSPsubj, MOT.pctMotion, ...
        conditions{iCond}, nPerms, seed);
end

%% Step 2b. Is clustering sensitive to a global inflation of connectivity?
fprintf('\n=== Step 2b: proportional-threshold invariance check ===\n');
IV = struct();
for iCond = 1:numel(conditions)
    IV.(conditions{iCond}) = rev01_threshold_invariance( ...
        dataDir, conditions{iCond}, 'HbO', [], seed);
end

%% Step 3. Pruning, motion and run exclusion counts
fprintf('\n=== Step 3: pruning report ===\n');
PR = rev01_pruning_report(dataDir);

%% Step 4. Participant characteristics and behavioural results
fprintf('\n=== Step 4: participants and behaviour ===\n');
BD = rev01_behavior_demographics(partDir, nPerms, seed);

%% Step 5. Save a compact result bundle
fprintf('\n=== Step 5: writing results ===\n');
R = struct();
R.generated    = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
R.nPerms       = nPerms;
R.seed         = seed;
[cL, cR] = get_channels_from_template('prefrontal');
R.sscIndices   = setdiff(1:24, [cL cR]);
R.longChannels = sort([cL cR]);

% Keep every field of the network metrics except the large per-threshold
% arrays. An allowlist was used here previously and silently dropped fields
% that the report then asked for, so the rule is now inverted: name only what
% is too big to keep.
dropNM = {'lambda', 'lambdaB', 'clust', 'clustW', 'SWI'};
fn = fieldnames(NM);
for k = 1:numel(fn)
    s = NM.(fn{k});
    sf = fieldnames(s);
    for j = 1:numel(sf)
        if ~ismember(sf{j}, dropNM)
            R.NM.(fn{k}).(sf{j}) = s.(sf{j});
        end
    end
end

% Same rule for the covariate analysis: keep everything except the model
% objects and the long-format table, which are large and not needed here.
dropQA = {'table', 'lmeUnadj', 'lmeAdj', 'anovaUnadj', 'anovaAdj', 'lrt', 'residClust'};
fn = fieldnames(QA);
for k = 1:numel(fn)
    s = QA.(fn{k});
    sf = fieldnames(s);
    for j = 1:numel(sf)
        if ~ismember(sf{j}, dropQA)
            R.QA.(fn{k}).(sf{j}) = s.(sf{j});
        end
    end
    R.QA.(fn{k}).anovaAdjTable   = dataset2cellSafe(s.anovaAdj);
    R.QA.(fn{k}).anovaUnadjTable = dataset2cellSafe(s.anovaUnadj);
    R.QA.(fn{k}).lrtTable        = dataset2cellSafe(s.lrt);
end

R.IV = IV;
R.PR = PR;
R.BD = BD;

save(fullfile(dataDir, 'rev01_results.mat'), '-struct', 'R');
rev01_write_report(R, fullfile(dataDir, 'rev01_results.txt'));

d = dir(fullfile(dataDir, 'rev01_results.mat'));
fprintf('\nSaved %s (%.0f kB) and rev01_results.txt\n', ...
        fullfile(dataDir, 'rev01_results.mat'), d.bytes/1024);
fprintf('Attach BOTH files to the conversation.\n');

%% ------------------------------------------------------------------------
function C = dataset2cellSafe(T)
% Convert an anova/compare result to a plain cell array so that it can be
% saved without the Statistics toolbox object wrapper.
try
    if isa(T, 'dataset'), T = dataset2table(T); end
    C = [T.Properties.VariableNames; table2cell(T)];
catch
    C = {};
end
end

% EOF
