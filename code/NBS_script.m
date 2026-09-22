%% NBS_script
clear; close all; format compact; clc

%% Create design matrix
%% NBS design matrix for 30 subjects x 3 time points
nSub = 30;
nTime = 3;

% Subject identity matrix
I = eye(nSub);

% Time coding:
% T1 = [1 0]
% T2 = [0 1]
% T3 = [0 0]   % reference level
T1 = repmat([1 0], nSub, 1);
T2 = repmat([0 1], nSub, 1);
T3 = repmat([0 0], nSub, 1);

% Design matrix: rows ordered as
% S1T1, S2T1, ..., S30T1,
% S1T2, S2T2, ..., S30T2,
% S1T3, S2T3, ..., S30T3
X = [
    I  T1
    I  T2
    I  T3
];

% Omnibus repeated-measures effect of time
contrast = [zeros(1,nSub) 1 1];   % F-test

% Exchange blocks: same subject ID across repeated measurements
exchange = repmat((1:nSub)', nTime, 1);

% Optional: package into GLM structure for NBS
GLM.X = X;
GLM.contrast = contrast;
GLM.exchange = exchange;
GLM.test = 'ftest';

% Check dimensions
disp(size(GLM.X))         % should be [90 32]
disp(size(GLM.contrast))  % should be [1 32]
disp(size(GLM.exchange))  % should be [90 1]

% renaming to NBS conventions
design = X;
save('..\data\design_NBS.mat', 'design', '-v7.3');
save('..\data\exchange_NBS.mat', 'exchange', '-v7.3');

%% post hoc paired contrasts
% T1 vs T2
c_T1_vs_T2 = [zeros(1,nSub)  1 -1];

% T1 vs T3
c_T1_vs_T3 = [zeros(1,nSub)  1  0];

% T2 vs T3
c_T2_vs_T3 = [zeros(1,nSub)  0  1];

% overall time effect: [zeros(1,30) 1 1] with F-test
% begin vs mid: [zeros(1,30) 1 -1]
% begin vs end: [zeros(1,30) 1 0]
% mid vs end: [zeros(1,30) 0 1]

%% Build 22x22x90 array for NBS from 3 time points (laparoscopic training)
% Assumptions:
% 1) Each .mat file contains a variable named zMatFDR
% 2) zMatFDR is a 30x1 cell array
% 3) Each cell contains one 24x24 connectivity matrix
% 4) Files are in the current folder, or provide full paths below

% File names in the desired stacking order:
% S1T1, S2T1, ..., S30T1, S1T2, ..., S30T2, S1T3, ..., S30T3
baseDir = '..\data\';
fileList = { ...
    fullfile(baseDir,'beginLapConnHbO.mat'), ...
    fullfile(baseDir,'midLapConnHbO.mat'), ...
    fullfile(baseDir,'endLapConnHbO.mat')};

nSub  = 30;
nTime = numel(fileList);
nObs  = nSub * nTime;

nodesToRemove = [3 14];
nodesAll      = 1:24;
nodesKeep     = setdiff(nodesAll, nodesToRemove);

nNodesFinal = numel(nodesKeep);   % should be 22

% Preallocate output array
connArray = nan(nNodesFinal, nNodesFinal, nObs);

idxObs = 0;

for idxTime = 1:nTime
    
    S = load(fileList{idxTime}, 'zMatFDR');
    
    if ~isfield(S, 'zMatFDR')
        error('File %s does not contain variable zMatFDR.', fileList{idxTime});
    end
    
    zMatFDR = S.zMatFDR;
    
    if ~iscell(zMatFDR) || numel(zMatFDR) ~= nSub
        error('In file %s, zMatFDR must be a %dx1 cell array.', fileList{idxTime}, nSub);
    end
    
    for idxSub = 1:nSub
        idxObs = idxObs + 1;
        
        thisMat = zMatFDR{idxSub};
        
        if ~ismatrix(thisMat) || ~isequal(size(thisMat), [24 24])
            error('Subject %d in file %s does not contain a 24x24 matrix.', ...
                idxSub, fileList{idxTime});
        end
        
        % Remove nodes 3 and 14
        % Short separation channels are [3, 14]
        % each fc matrix should be 22x22 long channels
        thisMat = thisMat(nodesKeep, nodesKeep);
        
        % Store in 3D array
        connArray(:, :, idxObs) = thisMat;
    end
end

% Verify final size
disp(size(connArray))   % should return [22 22 90]

% Renaming to NBS convention
Mat = connArray;

% Optional: save result
save('..\data\conn_Mat_HbO_NBS.mat', 'Mat', '-v7.3');

%% Build 22x22x90 array for NBS from 3 time points (resting-state)
% File names in the desired stacking order:
% S1T1, S2T1, ..., S30T1, S1T2, ..., S30T2, S1T3, ..., S30T3
baseDir = '..\data\';
fileList = { ...
    fullfile(baseDir,'beginRestingConnHbO.mat'), ...
    fullfile(baseDir,'midRestingConnHbO.mat'), ...
    fullfile(baseDir,'endRestingConnHbO.mat')};

nSub  = 30;
nTime = numel(fileList);
nObs  = nSub * nTime;

nodesToRemove = [3 14];
nodesAll      = 1:24;
nodesKeep     = setdiff(nodesAll, nodesToRemove);

nNodesFinal = numel(nodesKeep);   % should be 22

% Preallocate output array
connArray = nan(nNodesFinal, nNodesFinal, nObs);

idxObs = 0;

for idxTime = 1:nTime
    
    S = load(fileList{idxTime}, 'zMatFDR');
    
    if ~isfield(S, 'zMatFDR')
        error('File %s does not contain variable zMatFDR.', fileList{idxTime});
    end
    
    zMatFDR = S.zMatFDR;
    
    if ~iscell(zMatFDR) || numel(zMatFDR) ~= nSub
        error('In file %s, zMatFDR must be a %dx1 cell array.', fileList{idxTime}, nSub);
    end
    
    for idxSub = 1:nSub
        idxObs = idxObs + 1;
        
        thisMat = zMatFDR{idxSub};
        
        if ~ismatrix(thisMat) || ~isequal(size(thisMat), [24 24])
            error('Subject %d in file %s does not contain a 24x24 matrix.', ...
                idxSub, fileList{idxTime});
        end
        
        % Remove nodes 3 and 14
        % Short separation channels are [3, 14]
        % each fc matrix should be 22x22 long channels
        thisMat = thisMat(nodesKeep, nodesKeep);
        
        % Store in 3D array
        connArray(:, :, idxObs) = thisMat;
    end
end

% Verify final size
disp(size(connArray))   % should return [22 22 90]

% Renaming to NBS convention
Mat = connArray;

% Optional: save result
save('..\data\conn_Mat_HbO_NBS_resting.mat', 'Mat', '-v7.3');

% EOF