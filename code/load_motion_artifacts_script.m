% Script to load motion artifacts from all 6 groups, save percentage of
% time as motion artifacts

%% load_motion_artifacts
clear; close all; format compact; clc;
check_homer_path
saveFigs = false;

%% Directory List
dirDeriv{1} = '..\data\lap01\derivatives\homer';
dirDeriv{2} = '..\data\lap02\derivatives\homer';
dirDeriv{3} = '..\data\lap03\derivatives\homer';
dirDeriv{4} = '..\data\resting01\derivatives\homer';
dirDeriv{5} = '..\data\resting02\derivatives\homer';
dirDeriv{6} = '..\data\resting03\derivatives\homer';

% Load and rename
% load(fullfile(dirDeriv{6}, 'sub-29\nirs', 'sub-29_task-resting03_nirs.mat'))
% figure; stem(~output.misc.tIncAuto{1})

%% Load all subjects from all groups
nGroups = numel(dirDeriv);
nSubs   = 30;

tIncAutoAll = cell(nSubs, nGroups);          % stores output.misc.tIncAuto{1}
pctMotion   = nan(nSubs, nGroups);           % percentage of time with motion artifacts

for idxGroup = 1:nGroups
    
    % Extract dirName from path, e.g. '..\data\resting03\derivatives\homer' -> 'resting03'
    pathParts = strsplit(dirDeriv{idxGroup}, filesep);
    dirName = pathParts{3};   % {'..','data','resting03','derivatives','homer'}
    
    for idxSub = 1:nSubs
        
        subStr = sprintf('sub-%02d', idxSub);
        fileName = sprintf('%s_task-%s_nirs.mat', subStr, dirName);
        filePath = fullfile(dirDeriv{idxGroup}, [subStr '\nirs'], fileName);
        
        if ~isfile(filePath)
            warning('File not found: %s', filePath);
            continue
        end
        
        S = load(filePath, 'output');
        
        if ~isfield(S, 'output') || ...
           ~isstruct(S.output.misc) || ...
           ~isfield(S.output.misc, 'tIncAuto') || ...
           isempty(S.output.misc.tIncAuto) || ...
           numel(S.output.misc.tIncAuto) < 1
            warning('Missing output.misc.tIncAuto{1} in: %s', filePath);
            continue
        end
        
        % Save vector
        tIncAutoAll{idxSub, idxGroup} = S.output.misc.tIncAuto{1};
        
        % Percentage of time with motion artifacts
        pctMotion(idxSub, idxGroup) = 100 * sum(~S.output.misc.tIncAuto{1}) / numel(S.output.misc.tIncAuto{1});
    end
    fprintf('Group %d of %d done!\n', idxGroup, nGroups);
end

%% Plot with boxscatter
markerSize = 50;
lineWidth = 2;
figure('Color','w');
% training and resting sessions interleaved
hMotArt = boxscatter(pctMotion(:,[1 4 2 5 3 6]), 'ooxx^^', [0 0 1;0 0 1;0 1 0;0 1 0;1 0 0;1 0 0], markerSize, lineWidth); % blue o, green x, red ^
ylabel('Time with motion artifacts (%)');
xticks(1:nGroups);
xticklabels({'training','resting','training','resting','training','resting'});
box off
set(gca,'LineWidth',lineWidth,'TickDir','out', 'FontSize', 14);
set(hMotArt, 'Color','w','Units','inches', 'Position',[1 1 7.5 5.0]);
pngFile = '..\figures\motion_artifacts_per.png';
figFile = '..\figures\motion_artifacts_per.fig';
if saveFigs
    exportgraphics(hMotArt, pngFile, 'Resolution', 1200, 'BackgroundColor', 'white');
    savefig(hMotArt, figFile);
end

%% Optional save
save('..\data\motion_artifacts_summary.mat', 'tIncAutoAll', 'pctMotion', 'dirDeriv');

% EOF