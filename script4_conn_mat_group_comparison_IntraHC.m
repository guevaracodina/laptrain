%% Intrahemispherical correlations script 
% This is the fourth script to be run (4/4)
% First, you have to run conn_mat_group_level_script.m
% Then you need to run conn_mat_3_group_comparison.m
% Optional: run conn_mat_group_comparison_InterHC_script.m
%% Controls vs. Parkinson's Disease
clear; 
groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
disp('Beginning vs. Middle vs. End of shift Intrahemispheric correlation')

%% Interhemispherical connectivity HbO - Laparoscopic Training
group1Color = [1 0 0];      % Beginning HbO color
group2Color = [1 0 1];      % Middle HbO color
group3Color = [1 0.5 0];    % End HbO color
fileGroup1 = '..\data\beginLapConnHbO.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbO.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbO.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_lap_IntraHC_HbO';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbO
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

%% Interhemispherical connectivity HbR - Laparoscopic Training
clear;
groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
group1Color = [0 0 1];      % Beginning HbR color
group2Color = [0 1 1];      % Middle HbR color
group3Color = [0 0 0.5];    % End HbR color (Deep Navy)
fileGroup1 = '..\data\beginLapConnHbR.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbR.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbR.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_lap_IntraHC_HbR';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbR
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

%% Interhemispherical connectivity HbT - Laparoscopic Training
clear;
groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
group1Color = [0.00787401574803150,0.412505789717462,0.284113015284854];      % Beginning HbT color
group2Color = [0 1 0];      % Middle HbT color
group3Color = [0.5, 0.8, 0.2];   % End HbT color (chartreuse-olive)
fileGroup1 = '..\data\beginLapConnHbT.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midLapConnHbT.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endLapConnHbT.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_lap_IntraHC_HbT';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbR
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

%% Interhemispherical connectivity HbO - Resting-state
group1Color = [1 0 0];      % Beginning HbO color
group2Color = [1 0 1];      % Middle HbO color
group3Color = [1 0.5 0];    % End HbO color
fileGroup1 = '..\data\beginRestingConnHbO.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbO.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbO.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_resting_IntraHC_HbO';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbO
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

%% Interhemispherical connectivity HbR - Resting-state
clear;
groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
group1Color = [0 0 1];      % Beginning HbR color
group2Color = [0 1 1];      % Middle HbR color
group3Color = [0 0 0.5];    % End HbR color (Deep Navy)
fileGroup1 = '..\data\beginRestingConnHbR.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbR.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbR.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_resting_IntraHC_HbR';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbR
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

%% Interhemispherical connectivity HbT - Resting-state
clear;
groupNames = {'Begin' 'Mid' 'End'};
alpha = 0.05;
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
group1Color = [0.00787401574803150,0.412505789717462,0.284113015284854];      % Beginning HbT color
group2Color = [0 1 0];      % Middle HbT color
group3Color = [0.5, 0.8, 0.2];   % End HbT color (chartreuse-olive)
fileGroup1 = '..\data\beginRestingConnHbT.mat';     % CHANGE AS NEEDED!
fileGroup2 = '..\data\midRestingConnHbT.mat';       % CHANGE AS NEEDED!
fileGroup3 = '..\data\endRestingConnHbT.mat';       % CHANGE AS NEEDED!
saveFileName = 'conn_resting_IntraHC_HbT';          % CHANGE AS NEEDED!
% Perform inter-hemispherical connectivity (InterHC) for HbR
[zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha);

% EOF