%% Full connectivity matrix comparison
% Comparison between beginning, middle and end of shift (repeated measures ANOVA)
% This is the second script to be run (2/4)
% First, you have to run conn_mat_group_level_script.m
clear; close all; clc
groupNames = {'o00h' 'o12h' 'o24h'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
disp('Beginning vs. Middle vs. End of shift')

%% HbO - Laparoscopic training
fileGroup1 = '..\data\beginLapConnHbO.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbO.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbO.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_lap_GroupCompHbO';         % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

%% HbR - Laparoscopic training
clear; groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
fileGroup1 = '..\data\beginLapConnHbR.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbR.mat';          % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbR.mat';          % CHANGE AS NEEDED!
saveFileName = 'conn_lap_GroupCompHbR'; % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

%% HbT - Laparoscopic training
clear; groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
fileGroup1 = '..\data\beginLapConnHbT.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbT.mat';          % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbT.mat';          % CHANGE AS NEEDED!
saveFileName = 'conn_lap_GroupCompHbT'; % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

%% HbO - Resting-state
fileGroup1 = '..\data\beginRestingConnHbO.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbO.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbO.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_resting_GroupCompHbO';         % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

%% HbR - Resting-state
clear; groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
fileGroup1 = '..\data\beginRestingConnHbR.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbR.mat';          % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbR.mat';          % CHANGE AS NEEDED!
saveFileName = 'conn_resting_GroupCompHbR'; % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

%% HbT - Resting-state
clear; groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
fileGroup1 = '..\data\beginRestingConnHbT.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbT.mat';          % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbT.mat';          % CHANGE AS NEEDED!
saveFileName = 'conn_resting_GroupCompHbT'; % CHANGE AS NEEDED!
[zMat1Mean, zMat2Mean, zMat3Mean, ...
          P12, P13, P23, pFDR12, pFDR13, pFDR23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, true);

% EOF
