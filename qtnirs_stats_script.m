%% Summarizes the results of QT-NIRS 
clear; close all; format compact; clc;
saveFigs = false;
figFileLap01 = '..\data\qtnirs\lap01\QC_Group_good-channels_scanwise.fig';
figFileLap02 = '..\data\qtnirs\lap02\QC_Group_good-channels_scanwise.fig';
figFileLap03 = '..\data\qtnirs\lap03\QC_Group_good-channels_scanwise.fig';
channelCounts(:,1) = qtnirs_extract_good_ch_scanwise(figFileLap01);
channelCounts(:,3) = qtnirs_extract_good_ch_scanwise(figFileLap02);
channelCounts(:,5) = qtnirs_extract_good_ch_scanwise(figFileLap03);
figFileRest01 = '..\data\qtnirs\resting01\QC_Group_good-channels_scanwise.fig';
figFileRest02 = '..\data\qtnirs\resting02\QC_Group_good-channels_scanwise.fig';
figFileRest03 = '..\data\qtnirs\resting03\QC_Group_good-channels_scanwise.fig';
channelCounts(:,2) = qtnirs_extract_good_ch_scanwise(figFileRest01);
channelCounts(:,4) = qtnirs_extract_good_ch_scanwise(figFileRest02);
channelCounts(:,6) = qtnirs_extract_good_ch_scanwise(figFileRest03);
% disp('Group-level HQ scans (SCI>=0.60 and PSP>=0.08) above 70% of the recording time')
% disp('mean and std. dev.')
% disp(mean(channelCounts))
% disp(std(channelCounts))
hChan = figure('Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
imagesc(channelCounts, [0 30]); colorbar; colormap('gray');
title('Group-level HQ scans (SCI>=0.60 and PSP>=0.08) above 70% of the recording time')
set(gca, 'XTickLabel', {'training' 'resting' 'training' 'resting' 'training' 'resting'})
ylabel('Channel #')
set(gca, 'FontSize', 14)
pngFile = '..\figures\qtnirs_HQ_scans.png';
figFile = '..\figures\qtnirs_HQ_scans.fig';
if saveFigs
    exportgraphics(hChan, pngFile, 'Resolution', 1200, 'BackgroundColor', 'white');
    savefig(hChan, figFile);
end

%% Load Report Tables 
reportFileRest01 = '..\data\qtnirs\resting01\QC_reportTable_resting01.mat';
reportFileRest02 = '..\data\qtnirs\resting02\QC_reportTable_resting02.mat';
reportFileRest03 = '..\data\qtnirs\resting03\QC_reportTable_resting03.mat';
reportFileLap01 = '..\data\qtnirs\lap01\QC_reportTable_lap01.mat';
reportFileLap02 = '..\data\qtnirs\lap02\QC_reportTable_lap02.mat';
reportFileLap03 = '..\data\qtnirs\lap03\QC_reportTable_lap03.mat';
load(reportFileRest01)
reportTableRest01 = myReportTable; clear myReportTable;
load(reportFileRest02)
reportTableRest02 = myReportTable; clear myReportTable;
load(reportFileRest03)
reportTableRest03 = myReportTable; clear myReportTable;
load(reportFileLap01)
reportTableLap01 = myReportTable; clear myReportTable;
load(reportFileLap02)
reportTableLap02 = myReportTable; clear myReportTable;
load(reportFileLap03)
reportTableLap03 = myReportTable; clear myReportTable;

%% HQ scans
nScans = numel(reportTableLap01);
nChannels = size(channelCounts, 1);
markerSize = 50;
lineWidth = 2;
channelCountsNorm = 100*channelCounts/nScans;
hPer = figure('Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
hPer = boxscatter(channelCountsNorm, 'ooxx^^', [0 0 1;0 0 1;0 1 0;0 1 0;1 0 0;1 0 0], markerSize, lineWidth);
title('Percentage of channels with HQ scans (SCI>=0.60 and PSP>=0.08) above 70% of the recording time')
set(gca, 'XTickLabel', {'training' 'resting' 'training' 'resting' 'training' 'resting'})
ylabel('High-quality channels [%]')
set(gca, 'FontSize', 14)
set(hPer, 'Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
pngFile = '..\figures\qtnirs_HQ_scans_per.png';
figFile = '..\figures\qtnirs_HQ_scans_per.fig';
if saveFigs
    exportgraphics(hPer, pngFile, 'Resolution', 1200, 'BackgroundColor', 'white');
    savefig(hPer, figFile);
end

%% SCI
for iScans = 1:nScans
    SCIrest01(:,iScans) = mean(reportTableRest01(iScans).sci_array,2, 'omitnan');
    SCIrest02(:,iScans) = mean(reportTableRest02(iScans).sci_array,2, 'omitnan');
    SCIrest03(:,iScans) = mean(reportTableRest03(iScans).sci_array,2, 'omitnan');
    SCIlap01(:,iScans) = mean(reportTableLap01(iScans).sci_array,2, 'omitnan');
    SCIlap02(:,iScans) = mean(reportTableLap02(iScans).sci_array,2, 'omitnan');
    SCIlap03(:,iScans) = mean(reportTableLap03(iScans).sci_array,2, 'omitnan');
end

SCIall = [mean(SCIlap01,2), mean(SCIrest01,2), mean(SCIlap02,2), mean(SCIrest02,2),...
    mean(SCIlap03,2), mean(SCIrest03,2), ];

hSCI = figure('Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
hSCI = boxscatter(SCIall, 'ooxx^^', [0 0 1;0 0 1;0 1 0;0 1 0;1 0 0;1 0 0], markerSize, lineWidth);
title('Scalp-Coupling Index')
set(gca, 'XTickLabel', {'training' 'resting' 'training' 'resting' 'training' 'resting'})
ylabel('Scalp-Coupling Index')
set(gca, 'FontSize', 14)
set(hSCI, 'Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);

pngFile = '..\figures\qtnirs_SCI.png';
figFile = '..\figures\qtnirs_SCI.fig';
if saveFigs
    exportgraphics(hSCI, pngFile, 'Resolution', 1200, 'BackgroundColor', 'white');
    savefig(hSCI, figFile);
end

%% PSP 
for iScans = 1:nScans
    PSPrest01(:,iScans) = mean(reportTableRest01(iScans).power_array,2, 'omitnan');
    PSPrest02(:,iScans) = mean(reportTableRest02(iScans).power_array,2, 'omitnan');
    PSPrest03(:,iScans) = mean(reportTableRest03(iScans).power_array,2, 'omitnan');
    PSPlap01(:,iScans) = mean(reportTableLap01(iScans).power_array,2, 'omitnan');
    PSPlap02(:,iScans) = mean(reportTableLap02(iScans).power_array,2, 'omitnan');
    PSPlap03(:,iScans) = mean(reportTableLap03(iScans).power_array,2, 'omitnan');
end

PSPall = [mean(PSPlap01,2), mean(PSPrest01,2), mean(PSPlap02,2), mean(PSPrest02,2),...
    mean(PSPlap03,2), mean(PSPrest03,2), ];
hPSP = figure('Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
hPSP = boxscatter(PSPall, 'ooxx^^', [0 0 1;0 0 1;0 1 0;0 1 0;1 0 0;1 0 0], markerSize, lineWidth);
title('Peak Spectral Power')
set(gca, 'XTickLabel', {'training' 'resting' 'training' 'resting' 'training' 'resting'})
ylabel('Peak Spectral Power')
set(gca, 'FontSize', 14)
set(hPSP, 'Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);

pngFile = '..\figures\qtnirs_PSP.png';
figFile = '..\figures\qtnirs_PSP.fig';
if saveFigs
    exportgraphics(hPSP, pngFile, 'Resolution', 1200, 'BackgroundColor', 'white');
    savefig(hPSP, figFile);
end

disp('Quality assessment report done!')

% EOF