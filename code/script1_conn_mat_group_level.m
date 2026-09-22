%% Compute connectivity matrix for all groups (begin, middle and end of shift)
% This is the first script to be run (1/4)
% Works with SNIRF-compatible files, located in derivatives/homer Homer3 (v1.54.0)

%% channel by channel data for laparoscopic training
clear; close all; clc
check_homer_path
disp('Laparoscopic training computing connectivity matrices (channel by channel)')
beginDir = '..\data\lap01\derivatives\homer';
midDir = '..\data\lap02\derivatives\homer';
endDir = '..\data\lap03\derivatives\homer';

for iHb = 1:3               % 1=HbO, 2=HbR, 3=HbT
    % Beginning of shift Group
    [beginLap.rMatFDR, beginLap.zMatFDR, beginLap.pMatFDR, beginLap.pMat, beginLap.keepRun] = ...
        conn_mat_group_level(beginDir, iHb, 'beginLapConn');
    % End of shift Group
    [midLap.rMatFDR, midLap.zMatFDR, midLap.pMatFDR, midLap.pMat, midLap.keepRun] = ...
        conn_mat_group_level(midDir, iHb, 'midLapConn');
    % End of shift Group
    [endLap.rMatFDR, endLap.zMatFDR, endLap.pMatFDR, endLap.pMat, endLap.keepRun] = ...
        conn_mat_group_level(endDir, iHb, 'endLapConn');
    % Get Hb string label
    strContrast = get_Hb_string(iHb);
    % Save all fc data
    save(fullfile('..\data\', ['lap_all_significant_connections_' strContrast '.mat']));
end
disp('Laparoscopic training connectivity matrices (channel by channel) computed!')
clear; close all; clc

%% channel by channel data for resting state
disp('Resting-state  computing connectivity matrices (channel by channel)')
beginDir = '..\data\resting01\derivatives\homer';
midDir = '..\data\resting02\derivatives\homer';
endDir = '..\data\resting03\derivatives\homer';

for iHb = 1:3               % 1=HbO, 2=HbR, 3=HbT
    % Beginning of shift Group
    [beginResting.rMatFDR, beginResting.zMatFDR, beginResting.pMatFDR, beginResting.pMat, beginResting.keepRun] = ...
        conn_mat_group_level(beginDir, iHb, 'beginRestingConn');
    % End of shift Group
    [midResting.rMatFDR, midResting.zMatFDR, midResting.pMatFDR, midResting.pMat, midResting.keepRun] = ...
        conn_mat_group_level(midDir, iHb, 'midRestingConn');
    % End of shift Group
    [endResting.rMatFDR, endResting.zMatFDR, endResting.pMatFDR, endResting.pMat, endResting.keepRun] = ...
        conn_mat_group_level(endDir, iHb, 'endRestingConn');
    % Get Hb string label
    strContrast = get_Hb_string(iHb);
    % Save all fc data
    save(fullfile('..\data\', ['resting_all_significant_connections_' strContrast '.mat']));
end
disp('Resting-state connectivity matrices (channel by channel) computed!')

% EOF