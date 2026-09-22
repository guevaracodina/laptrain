function qc_summary_export(dataDir, outFile)
% QC_SUMMARY_EXPORT  Condense QT-NIRS report tables into a compact summary.
%
%   qc_summary_export()
%   qc_summary_export(dataDir)
%   qc_summary_export(dataDir, outFile)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Written to answer Associate Editor major concern 1d, which asks whether
%   the mid-shift clustering effect survives adjustment for data quality.
%   The QC_reportTable_*.mat files produced by QT-NIRS are about 142 MB each
%   and are therefore impractical to load repeatedly or to share. This
%   function reduces them to per-scan and per-channel summaries of the
%   scalp-coupling index (SCI) and the peak spectral power (PSP), together
%   with the scan-wise proportion of channels meeting the quality
%   criterion, and writes a single small .mat file.
%
%   INPUT
%       dataDir : char, optional
%           Path to the project 'data' folder. Default '..\data'.
%       outFile : char, optional
%           Output .mat file. Default fullfile(dataDir,'qc_summary.mat').
%
%   OUTPUT (variables saved in outFile)
%       SCIsubj      [nSub x 6] mean SCI per scan, long channels only
%       PSPsubj      [nSub x 6] mean PSP per scan, long channels only
%       SCIsubjAll   [nSub x 6] mean SCI per scan, all channels
%       PSPsubjAll   [nSub x 6] mean PSP per scan, all channels
%       SCIchan      [nCh x 6]  mean SCI per channel, averaged over scans
%       PSPchan      [nCh x 6]  mean PSP per channel, averaged over scans
%       goodChanPct  [nSub x 6] percentage of channels whose fraction of
%                               good 5 s windows reaches qualityThreshold
%       goodChanFrac {1 x 6}    [nCh x nSub] fraction of good windows
%       qualityThr   scalar     QT-NIRS master quality threshold
%       condNames    {1 x 6}    condition order of the columns
%       longChIdx    [1 x 22]   long-channel indices used
%
%   Column order is {'lap01','lap02','lap03','resting01','resting02','resting03'}.
%
%   Run once from the 'code' folder; afterwards only qc_summary.mat is needed.

if nargin < 1 || isempty(dataDir)
    dataDir = fullfile('..', 'data');
end
if nargin < 2 || isempty(outFile)
    outFile = fullfile(dataDir, 'qc_summary.mat');
end

condNames = {'lap01', 'lap02', 'lap03', 'resting01', 'resting02', 'resting03'};
nCond     = numel(condNames);

% Long channels of the 24-channel prefrontal Brite montage.
% Short-separation channels are indices 3 and 24; see get_channels_from_template.
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
longChIdx = sort([chIdxL, chIdxR]);

SCIsubj = []; PSPsubj = []; SCIsubjAll = []; PSPsubjAll = [];
SCIchan = []; PSPchan = []; goodChanPct = [];
goodChanFrac = cell(1, nCond);
qualityThr = NaN;

for iCond = 1:nCond

    reportFile = fullfile(dataDir, 'qtnirs', condNames{iCond}, ...
                          ['QC_reportTable_' condNames{iCond} '.mat']);
    if ~isfile(reportFile)
        error('qc_summary_export:missingFile', 'Not found: %s', reportFile);
    end

    fprintf('Loading %s ...\n', reportFile);
    S = load(reportFile, 'myReportTable');
    R = S.myReportTable;
    clear S

    nScans = numel(R);

    if iCond == 1
        nChan        = size(R(1).sci_array, 1);
        qualityThr   = R(1).thresholds.quality;
        SCIsubj      = nan(nScans, nCond);
        PSPsubj      = nan(nScans, nCond);
        SCIsubjAll   = nan(nScans, nCond);
        PSPsubjAll   = nan(nScans, nCond);
        SCIchan      = nan(nChan,  nCond);
        PSPchan      = nan(nChan,  nCond);
        goodChanPct  = nan(nScans, nCond);
    end

    sciChanTmp = nan(nChan, nScans);
    pspChanTmp = nan(nChan, nScans);
    fracTmp    = nan(nChan, nScans);

    for iScan = 1:nScans
        % sci_array and power_array are [nChannels x nWindows].
        % Average across the 5 s windows to obtain one value per channel.
        sciCh = mean(R(iScan).sci_array,   2, 'omitnan');
        pspCh = mean(R(iScan).power_array, 2, 'omitnan');

        sciChanTmp(:, iScan) = sciCh;
        pspChanTmp(:, iScan) = pspCh;

        % Scan-level summaries used as covariates in the revised analysis
        SCIsubjAll(iScan, iCond) = mean(sciCh, 'omitnan');
        PSPsubjAll(iScan, iCond) = mean(pspCh, 'omitnan');
        SCIsubj(iScan, iCond)    = mean(sciCh(longChIdx), 'omitnan');
        PSPsubj(iScan, iCond)    = mean(pspCh(longChIdx), 'omitnan');

        % good_combo_link(:,3) holds, per channel, the fraction of 5 s
        % windows in which both SCI and PSP exceed their thresholds.
        frac = R(iScan).good_combo_link(:, 3);
        fracTmp(1:numel(frac), iScan) = frac;
        goodChanPct(iScan, iCond) = 100 * mean(frac >= qualityThr);
    end

    SCIchan(:, iCond)   = mean(sciChanTmp, 2, 'omitnan');
    PSPchan(:, iCond)   = mean(pspChanTmp, 2, 'omitnan');
    goodChanFrac{iCond} = fracTmp;

    clear R
end

save(outFile, 'SCIsubj', 'PSPsubj', 'SCIsubjAll', 'PSPsubjAll', ...
              'SCIchan', 'PSPchan', 'goodChanPct', 'goodChanFrac', ...
              'qualityThr', 'condNames', 'longChIdx');

d = dir(outFile);
fprintf('Saved %s (%.1f kB)\n', outFile, d.bytes/1024);

end

% EOF
