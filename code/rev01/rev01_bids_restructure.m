function rev01_bids_restructure(srcRoot, dstRoot, partDir)
% REV01_BIDS_RESTRUCTURE  Reorganize the dataset to NIRS-BIDS (BIDS 1.8.0).
%
%   rev01_bids_restructure()
%   rev01_bids_restructure(srcRoot, dstRoot, partDir)
%
%   Neurophotonics NPH-260108-1, revision 01, Associate Editor major concern 2.
%
%   The released layout has six condition roots, each with sub-XX/nirs/, no
%   session level, no dataset_description.json and no sidecars. This function
%   writes a conformant tree in which the three visits are BIDS sessions:
%
%     <dstRoot>/dataset_description.json
%     <dstRoot>/participants.tsv, participants.json
%     <dstRoot>/README
%     <dstRoot>/sub-XX/ses-<00h|12h|24h>/nirs/
%         sub-XX_ses-<ses>_task-<lap|rest>_nirs.snirf
%         sub-XX_ses-<ses>_task-<lap|rest>_nirs.json
%         sub-XX_ses-<ses>_task-<lap|rest>_channels.tsv
%         sub-XX_ses-<ses>_task-<lap|rest>_optodes.tsv
%         sub-XX_ses-<ses>_task-<lap|rest>_coordsystem.json
%         sub-XX_ses-<ses>_task-lap_events.tsv        (training runs only)
%
%   Optode positions are read from the SNIRF probe rather than assumed, so the
%   _optodes.tsv reflects the acquisition geometry. Channel-wise MNI
%   coordinates and anatomical labels come from get_channel_info_prefrontal.
%
%   Step timings in residents_lap_train_tasks_timing.csv are converted to
%   per-run _events.tsv. Marker L is the video-synchronisation marker and is
%   written as onset 0 with trial_type 'video_start'; the remaining markers
%   become the 16 protocol steps. Cells holding '-' or left empty mean the
%   step was not reached and produce no event. Three cells in the source sheet
%   use a decimal point instead of a colon; these are repaired and logged.
%
%   AFTER RUNNING: validate the tree, e.g.
%       deno run -A jsr:@bids/validator <dstRoot>
%   and record the validator version and output in the manuscript.

if nargin < 1 || isempty(srcRoot), srcRoot = fullfile('..','data','BIDSsource'); end
if nargin < 2 || isempty(dstRoot), dstRoot = fullfile('..','data','BIDS_rev01'); end
if nargin < 3 || isempty(partDir), partDir = fullfile('..','data','participants'); end

conds = { 'lap01','ses-00h','lap'  ; 'lap02','ses-12h','lap'  ; 'lap03','ses-24h','lap'
          'resting01','ses-00h','rest'; 'resting02','ses-12h','rest'; 'resting03','ses-24h','rest' };
nSub = 30;

if ~isfolder(dstRoot), mkdir(dstRoot); end
info = get_channel_info_prefrontal();

%% dataset_description.json
dd = struct('Name', 'Resident physician mental workload assessed with fNIRS', ...
            'BIDSVersion', '1.8.0', ...
            'DatasetType', 'raw', ...
            'License', 'CC BY 4.0', ...
            'Authors', {{'Guevara, Edgar','Torres-Cuevas, Gerardo Enrique', ...
                         'Avalos Martinez, Javier','Kolosovas-Machuca, Eleazar Samuel', ...
                         'Martinez-Jimenez, Mario Aurelio'}}, ...
            'DatasetDOI', 'doi:10.5281/zenodo.15186569');
writejson(fullfile(dstRoot,'dataset_description.json'), dd);

fid = fopen(fullfile(dstRoot,'README'),'w');
fprintf(fid, ['Longitudinal fNIRS dataset of cognitive workload in laparoscopic training\n' ...
              'across a 24-hour clinical shift.\n\n' ...
              'Thirty fourth-year General Surgery residents performed a standardized\n' ...
              'laparoscopic training task (Origami Box Folding Exercise) and an\n' ...
              'eyes-closed resting-state recording at 0, 12 and 24 hours of a single\n' ...
              'continuous duty cycle. Sessions ses-00h, ses-12h and ses-24h correspond\n' ...
              'to those three visits. Tasks are task-lap and task-rest.\n\n' ...
              'The prefrontal montage has 24 channels; channels 3 (S03D01) and 24\n' ...
              '(S10D08) are short-separation channels.\n']);
fclose(fid);

%% participants.tsv
copyIfPresent(fullfile(partDir,'participants.tsv'), fullfile(dstRoot,'participants.tsv'));
pj = struct('age', struct('Description','Age of the participant','Units','years'), ...
            'sex', struct('Description','Biological sex', ...
                          'Levels', struct('m','male','f','female')), ...
            'handedness', struct('Description','Self-reported handedness', ...
                          'Levels', struct('r','right','l','left')), ...
            'PSCscore01', struct('Description','Procedure-Specific Checklist at 0 h, 0-19'), ...
            'PSCscore02', struct('Description','Procedure-Specific Checklist at 12 h, 0-19'), ...
            'PSCscore03', struct('Description','Procedure-Specific Checklist at 24 h, 0-19'), ...
            'GSRscore01', struct('Description', ...
               ['Global Rating Scale, 0-25. Recorded once per participant as a measure of ' ...
                'general operative skill; identical across the three visits by design.']));
writejson(fullfile(dstRoot,'participants.json'), pj);

%% Step timings, parsed once
[events, repairs] = readTimingEvents(fullfile(partDir,'residents_lap_train_tasks_timing.csv'), nSub);
if ~isempty(repairs)
    fprintf('Timing sheet: %d cell(s) repaired or treated as not reached.\n', numel(repairs));
end

%% Per run
nWritten = 0;
for iC = 1:size(conds,1)
    [oldTask, ses, task] = conds{iC,:};
    for iS = 1:nSub
        sub = sprintf('sub-%02d', iS);
        srcFile = fullfile(srcRoot, oldTask, sub, 'nirs', ...
                           sprintf('%s_task-%s_nirs.snirf', sub, oldTask));
        if ~isfile(srcFile)
            warning('rev01_bids_restructure:missing','Not found: %s', srcFile); continue
        end

        outDir = fullfile(dstRoot, sub, ses, 'nirs');
        if ~isfolder(outDir), mkdir(outDir); end
        base = sprintf('%s_%s_task-%s', sub, ses, task);

        copyfile(srcFile, fullfile(outDir, [base '_nirs.snirf']));

        snirf = SnirfLoad(srcFile);
        writeOptodes(fullfile(outDir,[base '_optodes.tsv']), snirf);
        writeCoordsystem(fullfile(outDir,[base '_coordsystem.json']));
        writeChannels(fullfile(outDir,[base '_channels.tsv']), snirf, info);
        writeNirsJson(fullfile(outDir,[base '_nirs.json']), snirf, task);

        if strcmp(task,'lap')
            iTime = find(strcmp(ses, {'ses-00h','ses-12h','ses-24h'}));
            writeEvents(fullfile(outDir,[base '_events.tsv']), events{iS, iTime});
        end
        nWritten = nWritten + 1;
    end
    fprintf('  %s -> %s / task-%s done\n', oldTask, ses, task);
end

fprintf(['\n%d run(s) written to %s\n' ...
         'Now run a BIDS validator on that folder and record its version and\n' ...
         'output in Section 3 of the manuscript.\n'], nWritten, dstRoot);
end

%% ------------------------------------------------------------------ helpers
function copyIfPresent(src, dst)
if isfile(src), copyfile(src, dst); else, warning('Not found: %s', src); end
end

function writejson(path, s)
fid = fopen(path, 'w'); assert(fid > 0, 'Cannot write %s', path);
fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true)); fclose(fid);
end

function writeOptodes(path, snirf)
% Optode positions come from the SNIRF probe, not from assumption.
p = snirf.probe;
if ~isempty(p.sourcePos3D), sp = p.sourcePos3D; else, sp = p.sourcePos2D; end
if ~isempty(p.detectorPos3D), dp = p.detectorPos3D; else, dp = p.detectorPos2D; end
fid = fopen(path,'w'); assert(fid>0,'Cannot write %s', path);
fprintf(fid, 'name\ttype\tx\ty\tz\n');
for i = 1:size(sp,1)
    fprintf(fid, 'S%d\tsource\t%.4f\t%.4f\t%.4f\n', i, sp(i,1), sp(i,2), col(sp,i,3));
end
for i = 1:size(dp,1)
    fprintf(fid, 'D%d\tdetector\t%.4f\t%.4f\t%.4f\n', i, dp(i,1), dp(i,2), col(dp,i,3));
end
fclose(fid);
end

function v = col(M, i, j)
if size(M,2) >= j, v = M(i,j); else, v = 0; end
end

function writeCoordsystem(path)
c = struct('NIRSCoordinateSystem', 'Other', ...
           'NIRSCoordinateUnits', 'mm', ...
           'NIRSCoordinateSystemDescription', ...
           ['Optode positions as stored in the SNIRF probe. Channel-wise MNI ' ...
            'coordinates and anatomical labels are given in the _channels.tsv ' ...
            'files and in optodes_MNI_coordinates.csv.']);
writejson(path, c);
end

function writeChannels(path, snirf, info)
ml = snirf.data(1).measurementList;
fid = fopen(path,'w'); assert(fid>0,'Cannot write %s', path);
fprintf(fid, 'name\ttype\tsource\tdetector\twavelength_nominal\tunits\tmni_x\tmni_y\tmni_z\tanatomical_label\n');
for k = 1:numel(ml)
    s = ml(k).sourceIndex; d = ml(k).detectorIndex;
    wl = snirf.probe.wavelengths(ml(k).wavelengthIndex);
    r = find(info.Source == s & info.Detector == d, 1);
    if isempty(r)
        x = NaN; y = NaN; z = NaN; lab = 'n/a';
    else
        x = info.X(r); y = info.Y(r); z = info.Z(r); lab = info.Label{r};
    end
    fprintf(fid, 'S%dD%d\tNIRSCWAMPLITUDE\tS%d\tD%d\t%g\tV\t%s\t%s\t%s\t%s\n', ...
            s, d, s, d, wl, num2str(x), num2str(y), num2str(z), lab);
end
fclose(fid);
end

function writeNirsJson(path, snirf, task)
t = snirf.data(1).time;
fs = 1 / mean(diff(t));
if strcmp(task,'lap')
    instr = 'Origami Box Folding Exercise performed in a laparoscopic box trainer, 5 min limit.';
else
    instr = 'Eyes-closed resting state, 6 min.';
end
j = struct('TaskName', task, ...
           'SamplingFrequency', fs, ...
           'NIRSChannelCount', numel(snirf.data(1).measurementList)/numel(snirf.probe.wavelengths), ...
           'NIRSSourceOptodeCount', size(snirf.probe.sourcePos2D,1), ...
           'NIRSDetectorOptodeCount', size(snirf.probe.detectorPos2D,1), ...
           'ShortChannelCount', 2, ...
           'Manufacturer', 'Artinis Medical Systems BV', ...
           'ManufacturersModelName', 'Brite MKII', ...
           'CapManufacturer', 'Artinis Medical Systems BV', ...
           'SoftwareFilters', 'n/a', ...
           'Instructions', instr);
writejson(path, j);
end

function writeEvents(path, ev)
fid = fopen(path,'w'); assert(fid>0,'Cannot write %s', path);
fprintf(fid, 'onset\tduration\ttrial_type\tstep_number\n');
if isempty(ev)
    fclose(fid); return
end
for k = 1:numel(ev.onset)
    fprintf(fid, '%.3f\tn/a\t%s\t%d\n', ev.onset(k), ev.label{k}, ev.step(k));
end
fclose(fid);
end

function [events, repairs] = readTimingEvents(fTim, nSub)
events  = cell(nSub, 3);
repairs = {};
if ~isfile(fTim)
    warning('rev01_bids_restructure:noTiming','Not found: %s', fTim); return
end
fid = fopen(fTim,'r','n','UTF-8'); c = onCleanup(@() fclose(fid));
lines = {}; while ~feof(fid), lines{end+1} = fgetl(fid); end %#ok<AGROW>

hdr = splitCSV(lines{1});
stepNames = hdr(2:end);

for iL = 3:numel(lines)
    L = lines{iL};
    if isempty(L) || ~ischar(L), continue; end
    f = splitCSV(L);
    tok = regexp(f{1}, 'sub\s*(\d+)\s*task\s*lap(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    iS = str2double(tok{1}); iT = str2double(tok{2});
    if isnan(iS) || isnan(iT) || iS > nSub || iT > 3, continue; end

    onset = []; label = {}; step = [];
    for j = 2:min(numel(f), numel(hdr))
        t = strtrim(f{j});
        if isempty(t) || strcmp(t,'-'), continue; end
        m = regexp(t, '^(\d+):(\d{2})$', 'tokens', 'once');
        if isempty(m)
            m = regexp(t, '^(\d+)\.(\d{2})$', 'tokens', 'once');
            if ~isempty(m)
                repairs{end+1} = sprintf('%s step %d: ''%s'' read as %s:%s', ...
                                         strtrim(f{1}), j-1, t, m{1}, m{2}); %#ok<AGROW>
            end
        end
        if isempty(m), continue; end
        sec = str2double(m{1})*60 + str2double(m{2});
        nm  = regexprep(stepNames{j-1}, '^\d+\s*', '');
        nm  = regexprep(lower(strtrim(nm)), '[^a-z0-9]+', '_');
        onset(end+1) = sec;                                   %#ok<AGROW>
        step(end+1)  = j-1;                                   %#ok<AGROW>
        if j == 2, label{end+1} = 'video_start';              %#ok<AGROW>
        else,      label{end+1} = nm; end                     %#ok<AGROW>
    end
    events{iS, iT} = struct('onset', onset, 'label', {label}, 'step', step);
end
end

function f = splitCSV(L)
f = {}; cur = ''; inQ = false;
for k = 1:numel(L)
    ch = L(k);
    if ch == '"', inQ = ~inQ;
    elseif ch == ',' && ~inQ, f{end+1} = cur; cur = ''; %#ok<AGROW>
    else, cur(end+1) = ch; %#ok<AGROW>
    end
end
f{end+1} = cur;
end

% EOF
