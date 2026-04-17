%% qtnirs: quality assessment of raw fNIRS data
clear all; close all; format compact; clc
check_homer_path
qtFolder= 'C:\Edgar\Dropbox\CIACYT\Projects\2023\fNIRS_mental_workload\code\qt_nirs';
addpath(qtFolder)
currFolder = pwd;
cd(qtFolder)
setpaths
cd(currFolder)

%% QT-NIRS parameters
bpFmin = 0.5; bpFmax = 2.5;
windowSec = 5;
windowOverlap = 0;
quality_threshold = 0.7;
sciThreshold = 0.6;
pspThreshold = 0.08;

%% Run QT-NIRS
qtnirs('..\data\qtnirs\lap01',...
                'freqCut',[bpFmin, bpFmax],...
                'window',windowSec,...
                'overlap',windowOverlap,...
                'qualityThreshold',quality_threshold,...
                'sciThreshold', sciThreshold,...
                'pspThreshold', pspThreshold,...
                'conditionsMask','all',...
                'dodFlag',0,...
                'guiFlag',0);
cd(currFolder)