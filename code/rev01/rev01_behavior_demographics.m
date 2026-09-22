function B = rev01_behavior_demographics(partDir, nPerms, seed)
% REV01_BEHAVIOR_DEMOGRAPHICS  Participant table and behavioural results.
%
%   B = rev01_behavior_demographics(partDir, nPerms, seed)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Addresses Associate Editor minor concerns 1 and 3.
%
%   Reads participants.tsv (participant_id, age, sex, handedness,
%   PSCscore01..03, GSRscore01..03) and residents_lap_train_tasks_timing.csv,
%   and derives:
%     - the participant characteristics table requested in concern 3;
%     - PSC at each time point, with repeated-measures permutation tests;
%     - task completion, defined as the number of the 17 protocol steps that
%       carry a timestamp, and the time of the last completed step.
%
%   IMPORTANT, see the warnings this function prints:
%     * GRS (column GSRscore01..03) is IDENTICAL at all three time points for
%       every subject, in participants.tsv and in the source workbook. It was
%       recorded once per subject, not per session, so it cannot be reported
%       as a repeated measure. The function reports it as a single
%       between-subject score and flags the fact.
%     * The timing sheet mixes '-' and empty cells for steps not reached, and
%       three cells use a decimal point instead of a colon (2.41, 4.56, 2.55).
%       Both are handled here, and every repair is listed.
%
%   INPUT
%       partDir  folder with participants.tsv and the timing CSV.
%                Default '..\data\participants'.
%       nPerms   permutations (default 5000)
%       seed     rng seed (default 42)

if nargin < 1 || isempty(partDir), partDir = fullfile('..', 'data', 'participants'); end
if nargin < 2 || isempty(nPerms),  nPerms = 5000; end
if nargin < 3 || isempty(seed),    seed   = 42;   end
rng(seed);

B = struct();

%% Participants
fPart = fullfile(partDir, 'participants.tsv');
assert(isfile(fPart), 'Not found: %s', fPart);
T = readtable(fPart, 'FileType', 'text', 'Delimiter', '\t', ...
              'VariableNamingRule', 'preserve');
T = T(~cellfun(@isempty, strtrim(string(T.participant_id))), :);
B.participants = T;
B.n = height(T);

age = double(T.age);
B.ageMean = mean(age, 'omitnan');
B.ageSD   = std(age,  'omitnan');
B.ageMin  = min(age);
B.ageMax  = max(age);

sex  = categorical(strtrim(string(T.sex)));
hand = categorical(strtrim(string(T.handedness)));
B.sexCategories  = categories(sex);   B.sexCounts  = countcats(sex);
B.handCategories = categories(hand);  B.handCounts = countcats(hand);

fprintf('\nParticipants: N = %d\n', B.n);
fprintf('  Age %.1f (%.1f) years, range %d to %d\n', B.ageMean, B.ageSD, B.ageMin, B.ageMax);
fprintf('  Sex: %s\n', joinCounts(B.sexCategories, B.sexCounts));
fprintf('  Handedness: %s\n', joinCounts(B.handCategories, B.handCounts));

% Variables requested by the Associate Editor that are absent from the data
missingVars = {};
for v = {'specialty', 'cap size', 'prior laparoscopic experience'}
    missingVars{end+1} = v{1}; %#ok<AGROW>
end
B.missingDemographics = missingVars;
warning('rev01_behavior_demographics:missingVars', ...
        ['participants.tsv holds no %s. These were requested in minor ' ...
         'concern 3 and must either be supplied or declared unavailable.'], ...
        strjoin(missingVars, ', '));

%% PSC across time points
PSC = [double(T.PSCscore01), double(T.PSCscore02), double(T.PSCscore03)];
B.PSCData = PSC;
B.PSCMean = mean(PSC, 1, 'omitnan');
B.PSCSD   = std(PSC, 0, 1, 'omitnan');
[B.PSCP, B.PSCDiff, B.PSCF] = permutationTestRM3(PSC, nPerms);
B.PSCP12 = permutationTestRMpair(PSC(:,1), PSC(:,2), nPerms);
B.PSCP13 = permutationTestRMpair(PSC(:,1), PSC(:,3), nPerms);
B.PSCP23 = permutationTestRMpair(PSC(:,2), PSC(:,3), nPerms);

fprintf('\nPSC  0 h / 12 h / 24 h = %.2f (%.2f) / %.2f (%.2f) / %.2f (%.2f), global p = %.4f\n', ...
        B.PSCMean(1), B.PSCSD(1), B.PSCMean(2), B.PSCSD(2), B.PSCMean(3), B.PSCSD(3), B.PSCP);

%% GRS: verify whether it actually varies across time
GRS = [double(T.GSRscore01), double(T.GSRscore02), double(T.GSRscore03)];
B.GRSData = GRS;
B.GRSConstantAcrossTime = all(all(GRS == GRS(:,1), 2));

if B.GRSConstantAcrossTime
    B.GRSScore = GRS(:,1);
    B.GRSMean  = mean(B.GRSScore, 'omitnan');
    B.GRSSD    = std(B.GRSScore,  'omitnan');
    B.GRSMin   = min(B.GRSScore);
    B.GRSMax   = max(B.GRSScore);
    warning('rev01_behavior_demographics:grsConstant', ...
            ['GRS is identical at all three time points for every subject. ' ...
             'It is a single between-subject score and cannot be tested ' ...
             'across time. Reported as one value per participant.']);
    fprintf('GRS  single score per subject: %.2f (%.2f), range %d to %d\n', ...
            B.GRSMean, B.GRSSD, B.GRSMin, B.GRSMax);
else
    B.GRSMean = mean(GRS, 1, 'omitnan');
    B.GRSSD   = std(GRS, 0, 1, 'omitnan');
    [B.GRSP, B.GRSDiff, B.GRSF] = permutationTestRM3(GRS, nPerms);
    B.GRSP12 = permutationTestRMpair(GRS(:,1), GRS(:,2), nPerms);
    B.GRSP13 = permutationTestRMpair(GRS(:,1), GRS(:,3), nPerms);
    B.GRSP23 = permutationTestRMpair(GRS(:,2), GRS(:,3), nPerms);
end

%% Task completion from the step timing sheet
fTim = fullfile(partDir, 'residents_lap_train_tasks_timing.csv');
if isfile(fTim)
    [B.stepsCompleted, B.lastStepSec, B.timingRepairs] = readTiming(fTim, B.n);

    B.completionMean = mean(B.stepsCompleted, 1, 'omitnan');
    B.completionSD   = std(B.stepsCompleted, 0, 1, 'omitnan');
    [B.completionP, B.completionDiff, B.completionF] = ...
        permutationTestRM3(B.stepsCompleted, nPerms);
    B.completionP12 = permutationTestRMpair(B.stepsCompleted(:,1), B.stepsCompleted(:,2), nPerms);
    B.completionP13 = permutationTestRMpair(B.stepsCompleted(:,1), B.stepsCompleted(:,3), nPerms);
    B.completionP23 = permutationTestRMpair(B.stepsCompleted(:,2), B.stepsCompleted(:,3), nPerms);

    B.lastStepMean = mean(B.lastStepSec, 1, 'omitnan');
    B.lastStepSD   = std(B.lastStepSec, 0, 1, 'omitnan');
    B.lastStepP    = permutationTestRM3(B.lastStepSec, nPerms);

    fprintf('Completion (steps of 17) 0 h / 12 h / 24 h = %.2f / %.2f / %.2f, global p = %.4f\n', ...
            B.completionMean, B.completionP);
    fprintf('Time of last completed step (s) = %.0f / %.0f / %.0f, global p = %.4f\n', ...
            B.lastStepMean, B.lastStepP);
    if ~isempty(B.timingRepairs)
        fprintf('Timing cells repaired or treated as not reached: %d (see B.timingRepairs)\n', ...
                numel(B.timingRepairs));
    end
else
    warning('rev01_behavior_demographics:noTiming', 'Not found: %s', fTim);
end

end

%% ------------------------------------------------------------------------
function s = joinCounts(cats, cnts)
parts = arrayfun(@(i) sprintf('%s = %d (%.1f%%)', cats{i}, cnts(i), ...
                              100*cnts(i)/sum(cnts)), 1:numel(cats), ...
                 'UniformOutput', false);
s = strjoin(parts, ', ');
end

%% ------------------------------------------------------------------------
function [steps, lastSec, repairs] = readTiming(fTim, nSub)
% Parse the step timing sheet into [nSub x 3] completion counts.
% Row 1 holds the step names, row 2 the SNIRF markers, rows 3+ the runs.
% Column 2 is the 'Start' marker L and is excluded from the step count.

fid = fopen(fTim, 'r', 'n', 'UTF-8');
c = onCleanup(@() fclose(fid));
lines = {};
while ~feof(fid), lines{end+1} = fgetl(fid); end %#ok<AGROW>

steps   = nan(nSub, 3);
lastSec = nan(nSub, 3);
repairs = {};

for iLine = 3:numel(lines)
    L = lines{iLine};
    if isempty(L) || ~ischar(L), continue; end
    f = splitCSV(L);
    if isempty(f) || isempty(strtrim(f{1})), continue; end

    tok = regexp(f{1}, 'sub\s*(\d+)\s*task\s*lap(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    iSub  = str2double(tok{1});
    iTime = str2double(tok{2});
    if isnan(iSub) || isnan(iTime) || iSub > nSub || iTime > 3, continue; end

    sec = nan(1, 17);
    for j = 2:min(numel(f), 18)
        t = strtrim(f{j});
        if isempty(t) || strcmp(t, '-')
            if isempty(t)
                repairs{end+1} = sprintf('%s step %d: empty cell treated as not reached', ...
                                         strtrim(f{1}), j-1); %#ok<AGROW>
            end
            continue
        end
        m = regexp(t, '^(\d+):(\d{2})$', 'tokens', 'once');
        if isempty(m)
            % Decimal point used instead of a colon, e.g. 2.41
            m = regexp(t, '^(\d+)\.(\d{2})$', 'tokens', 'once');
            if ~isempty(m)
                repairs{end+1} = sprintf('%s step %d: ''%s'' read as %s:%s', ...
                                         strtrim(f{1}), j-1, t, m{1}, m{2}); %#ok<AGROW>
            end
        end
        if ~isempty(m)
            sec(j-1) = str2double(m{1})*60 + str2double(m{2});
        else
            repairs{end+1} = sprintf('%s step %d: unparsed ''%s''', ...
                                     strtrim(f{1}), j-1, t); %#ok<AGROW>
        end
    end

    done = ~isnan(sec(2:end));              % step 1 is the start marker
    steps(iSub, iTime)   = sum(done);
    if any(done)
        lastSec(iSub, iTime) = max(sec(2:end));   % max ignores NaN
    end
end
end

%% ------------------------------------------------------------------------
function f = splitCSV(L)
% Split one CSV line, honouring double-quoted fields.
f = {}; cur = ''; inQ = false;
for k = 1:numel(L)
    ch = L(k);
    if ch == '"'
        inQ = ~inQ;
    elseif ch == ',' && ~inQ
        f{end+1} = cur; cur = ''; %#ok<AGROW>
    else
        cur(end+1) = ch; %#ok<AGROW>
    end
end
f{end+1} = cur;
end

% EOF
