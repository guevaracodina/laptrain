function P = rev01_pruning_report(dataDir, saveFile, snirfRoot)
% REV01_PRUNING_REPORT  Channel pruning, motion and run exclusion counts.
%
%   P = rev01_pruning_report(dataDir)
%   P = rev01_pruning_report(dataDir, saveFile)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Answers Reviewer 1: how many channels were pruned, how many trials
%   rejected, and was any subject or block excluded outright.
%
%   REVISION NOTE. An earlier version of this function counted pruning across
%   every measurement in the raw SNIRF file. The Brite MKII records 96
%   source-detector combinations, of which only the 24 of the prefrontal
%   montage enter output.dc and the connectivity analysis, so that count was
%   reported against the wrong denominator and could exceed 24. Pruning is now
%   counted only among the 24 montage channels, matched by source-detector
%   pair, and both numbers are returned so the difference is visible.
%
%   INPUT
%       dataDir    project 'data' folder. Expects
%                  <dataDir>\<cond>\derivatives\homer\sub-XX\nirs\*.mat
%       saveFile   optional .mat path for the result
%       snirfRoot  folder holding the raw SNIRF files, used to recover the raw
%                  measurement list. Default <dataDir>\BIDSsource.
%
%   WHY THE SNIRF IS NEEDED. The saved Homer3 object is a ProcResultClass and
%   carries the processed dc but not the raw data, so out.data does not exist.
%   mlActAuto is indexed on the RAW measurement list (96 channels x 2
%   wavelengths on this system), which only the SNIRF file holds. The raw list
%   is read once per condition, not once per run.
%
%   OUTPUT struct P with [nSub x 6] matrices, columns ordered
%       {lap01 lap02 lap03 resting01 resting02 resting03}
%       nChanPruned     pruned among the 24 montage channels
%       nChanPrunedRaw  pruned among all raw channels
%       nChanMontage    montage channels found (should be 24)
%       nChanRaw        raw channels present in the file
%       pctMotion, retainedSec, keptRun, missingFile

REPORT_VERSION = '2026-09-18e';
fprintf('rev01_pruning_report version %s\n', REPORT_VERSION);

if nargin < 1 || isempty(dataDir), dataDir = fullfile('..', 'data'); end
if nargin < 3 || isempty(snirfRoot), snirfRoot = fullfile(dataDir, 'BIDSsource'); end

condNames = {'lap01', 'lap02', 'lap03', 'resting01', 'resting02', 'resting03'};
nCond = numel(condNames);
nSub  = 30;
minRecordingTime = 120;      % seconds, as in conn_mat_subject_level.m

% The 24 montage channels, as source-detector pairs, in Homer3 order
info = get_channel_info_prefrontal();
montage = [info.Source, info.Detector];
nMont   = height(info);        % 24

P.version        = REPORT_VERSION;
P.condNames      = condNames;
P.nChanPruned    = nan(nSub, nCond);
P.nChanPrunedRaw = nan(nSub, nCond);
P.nChanMatched   = nan(nSub, nCond);
P.nChanMontage   = nan(nSub, nCond);
P.nChanRaw       = nan(nSub, nCond);
P.pctMotion      = nan(nSub, nCond);
P.retainedSec    = nan(nSub, nCond);
P.keptRun        = false(nSub, nCond);
P.missingFile    = false(nSub, nCond);

warnedLen = false;
reportedShape = false;
for iCond = 1:nCond

    % ---- raw measurement list for this condition, read once from the SNIRF --
    srcAll = []; detAll = []; pairs = []; nRaw = 0; inMontage = [];
    snirfFile = findSnirf(snirfRoot, condNames{iCond});
    if isempty(snirfFile)
        fprintf('  %s: no SNIRF found under %s (informational only)\n', ...
                condNames{iCond}, snirfRoot);
    else
        sn  = SnirfLoad(snirfFile);
        ml  = sn.data(1).measurementList;
        srcAll = [ml.sourceIndex]';
        detAll = [ml.detectorIndex]';
        pairs  = unique([srcAll detAll], 'rows', 'stable');
        nRaw   = size(pairs, 1);
        inMontage = ismember(pairs, montage, 'rows');
        fprintf('  %s: raw list from %s (%d measurements, %d channels, %d in montage)\n', ...
                condNames{iCond}, snirfFile, numel(ml), nRaw, sum(inMontage));
    end

    for iSub = 1:nSub
        subStr   = sprintf('sub-%02d', iSub);
        fileName = sprintf('%s_task-%s_nirs.mat', subStr, condNames{iCond});
        filePath = fullfile(dataDir, condNames{iCond}, 'derivatives', 'homer', ...
                            subStr, 'nirs', fileName);

        if ~isfile(filePath)
            P.missingFile(iSub, iCond) = true;
            continue
        end

        S = load(filePath, 'output');
        out = S.output;

        % --- channel pruning, restricted to the montage ----------------------
        % mlActAuto{1} is [nMeas x 4] in Homer's measurement-list convention:
        %   column 1 source index
        %   column 2 detector index
        %   column 3 ACTIVE FLAG, 0 = pruned
        %   column 4 wavelength index
        % A channel counts as pruned if it was deactivated at either
        % wavelength. Rows are matched to the montage by source-detector pair
        % rather than by position, so the count cannot be thrown off by a
        % different measurement ordering.
        if isfield(out.misc, 'mlActAuto') && ~isempty(out.misc.mlActAuto)
            a = out.misc.mlActAuto{1};
            if ~reportedShape
                fprintf(['    mlActAuto{1} is %s; columns 3 and 4 take values %s ' ...
                         'and %s; dc list has %d entries\n'], mat2str(size(a)), ...
                        mat2str(unique(a(:,3))'), mat2str(unique(a(:,4))'), ...
                        numel(out.dc.measurementList));
                reportedShape = true;
            end
            if size(a,2) >= 4
                nPruned = 0; nSeen = 0;
                for m = 1:nMont
                    sel = a(:,1) == montage(m,1) & a(:,2) == montage(m,2);
                    if ~any(sel), continue; end
                    nSeen = nSeen + 1;
                    if ~all(a(sel,3)), nPruned = nPruned + 1; end
                end
                P.nChanPruned(iSub, iCond)  = nPruned;
                P.nChanMatched(iSub, iCond) = nSeen;
                P.nChanRaw(iSub, iCond)     = size(unique(a(:,1:2),'rows'),1);
            elseif ~warnedLen
                warning('rev01_pruning_report:mlShape', ...
                        'mlActAuto{1} is %s, not the expected [nMeas x 4].', mat2str(size(a)));
                warnedLen = true;
            end
        end

        % --- motion artifacts and retained time ------------------------------
        if isfield(out.misc, 'tIncAuto') && ~isempty(out.misc.tIncAuto)
            tInc = out.misc.tIncAuto{1};
            P.pctMotion(iSub, iCond) = 100 * sum(~tInc) / numel(tInc);
            t  = out.dc.GetTime();
            ts = mean(diff(t));
            P.retainedSec(iSub, iCond) = sum(tInc) * ts;
            P.keptRun(iSub, iCond) = P.retainedSec(iSub, iCond) >= minRecordingTime;
        end
    end
    fprintf('  pruning report: %s done\n', condNames{iCond});
end

%% Summary
P.meanChanPruned    = mean(P.nChanPruned, 1, 'omitnan');
P.sdChanPruned      = std(P.nChanPruned, 0, 1, 'omitnan');
P.maxChanPruned     = max(P.nChanPruned, [], 1);
P.meanChanPrunedRaw = mean(P.nChanPrunedRaw, 1, 'omitnan');
P.nRunsExcluded     = sum(~P.keptRun & ~P.missingFile, 1);
P.meanRetainedSec   = mean(P.retainedSec, 1, 'omitnan');
P.minRecordingTime  = minRecordingTime;
P.subjectsFullyExcluded = find(all(~P.keptRun, 2))';

matched = unique(P.nChanMatched(~isnan(P.nChanMatched)));
fprintf('\nMontage channels found in mlActAuto per run: %s of %d expected\n', ...
        mat2str(matched(:)'), nMont);
if ~isequal(matched(:)', nMont)
    warning('rev01_pruning_report:matchCount', ...
            'Not every montage channel was matched; counts may understate pruning.');
end

fprintf('\nChannels pruned per scan, of the %d montage channels:\n', nMont);
fprintf('  %-10s %-24s %-8s %s\n', 'condition', 'pruned, mean (SD)', 'max', 'runs below min');
for iCond = 1:nCond
    fprintf('  %-10s %8.2f (%5.2f)%9s %-8d %d\n', condNames{iCond}, ...
        P.meanChanPruned(iCond), P.sdChanPruned(iCond), '', ...
        P.maxChanPruned(iCond), P.nRunsExcluded(iCond));
end

fprintf('\nMean retained recording time per condition (s): %s\n', ...
        num2str(P.meanRetainedSec, '%.1f  '));
if isempty(P.subjectsFullyExcluded)
    fprintf('No subject was excluded in every condition.\n');
else
    fprintf('Subjects excluded in every condition: %s\n', num2str(P.subjectsFullyExcluded));
end

if nargin > 1 && ~isempty(saveFile)
    save(saveFile, '-struct', 'P');
    fprintf('Saved %s\n', saveFile);
end

end

%% ------------------------------------------------------------------------
function f = findSnirf(snirfRoot, condName)
% Locate any SNIRF for this condition; the montage is identical across runs,
% so the first one found is sufficient to recover the raw measurement list.
f = '';
cands = { fullfile(snirfRoot, condName, 'sub-01', 'nirs', ...
                   sprintf('sub-01_task-%s_nirs.snirf', condName)) };
for k = 1:numel(cands)
    if isfile(cands{k}), f = cands{k}; return; end
end
d = dir(fullfile(snirfRoot, condName, '**', '*.snirf'));
if ~isempty(d), f = fullfile(d(1).folder, d(1).name); return; end
d = dir(fullfile(snirfRoot, '**', sprintf('*task-%s_nirs.snirf', condName)));
if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
end

% EOF
