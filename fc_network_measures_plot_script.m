%% Load data
clear; close all; clc;
addpath(genpath('C:\Edgar\Dropbox\Matlab\BCT'))
addpath('C:\Edgar\Dropbox\Matlab\circularGraph')
addpath('C:\Edgar\Dropbox\Matlab\shadedErrorBar')
% Load data
load('..\data\lap_all_significant_connections_HbT.mat');
rng(42); % For repeatability
printFigures = false;
symmColorZero = true;

%% Prepare data
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];
% meanMat22 is your 22×22 mean connectivity matrix.
meanMatBegLap = meanMatFromCell(beginLap.zMatFDR, channelIdx);
meanMatMidLap = meanMatFromCell(midLap.zMatFDR, channelIdx);
meanMatEndLap = meanMatFromCell(endLap.zMatFDR, channelIdx);
channelLabels = {'Frontal Inf Tri R'
'Frontal Inf Tri R'
'SSC'
'Frontal Inf Tri R'
'Frontal Sup R'
'Frontal Mid R'
'Frontal Mid R'
'Frontal Sup R'
'Frontal Sup R'
'Frontal Sup R'
'Frontal Sup R'
'Frontal Sup Medial L'
'Frontal Sup Medial L'
'Frontal Sup R'
'Frontal Sup L'
'Frontal Mid L'
'Frontal Sup L'
'Frontal Mid L'
'Frontal Inf Tri L'
'Frontal Sup Mid L'
'Frontal Mid L'
'Frontal Inf Tri L'
'Frontal Inf Tri L'
'SSC'
};
channelLabels = channelLabels(channelIdx');

%% Plot circular graph of mean connectivity matrices
mn = min([meanMatBegLap(:); meanMatMidLap(:); meanMatEndLap(:)]); 
mx = max([meanMatBegLap(:); meanMatMidLap(:); meanMatEndLap(:)]);

%% Color limits for means (row1)
if symmColorZero
    limMean = max(abs([mn, mx]));
    mn = -limMean;
    mx =  limMean;
end

hFig1a = figure; set(hFig1a, 'color', 'w', 'name', 'Begin')
circularGraph(meanMatBegLap, 'Colormap',ioi_get_colormap('redbluecmap', 22),'MinValue',mn, 'MaxValue',mx,'Label',channelLabels);
hFig1b = figure; set(hFig1b, 'color', 'w', 'name', 'Middle')
circularGraph(meanMatMidLap, 'Colormap',ioi_get_colormap('redbluecmap', 22),'MinValue',mn, 'MaxValue',mx,'Label',channelLabels);
hFig1c = figure; set(hFig1c, 'color', 'w', 'name', 'End')
circularGraph(meanMatEndLap, 'Colormap',ioi_get_colormap('redbluecmap', 22),'MinValue',mn, 'MaxValue',mx,'Label',channelLabels);

% hFig1d = figure; set(hFig1c, 'color', 'w', 'name', 'significant')
% circularGraph(-log(pFDR), 'Colormap',ioi_get_colormap('redbluecmap', 22),'MinValue',mn, 'MaxValue',mx,'Label',channelLabels);

% Specify window units
set(hFig1a, 'units', 'inches')
% Change figure and paper size
set(hFig1a, 'Position', [0.1 0.1 10 6])
set(hFig1a, 'PaperPosition', [0.1 0.1 10 6])

% Specify window units
set(hFig1b, 'units', 'inches')
% Change figure and paper size
set(hFig1b, 'Position', [0.1 0.1 10 6])
set(hFig1b, 'PaperPosition', [0.1 0.1 10 6])

% Specify window units
set(hFig1c, 'units', 'inches')
% Change figure and paper size
set(hFig1c, 'Position', [0.1 0.1 10 6])
set(hFig1c, 'PaperPosition', [0.1 0.1 10 6])

if printFigures
    % Save as PNG
    print(hFig1a, '-dpng', fullfile('..\figures\', 'mean_conn_mat_lap_beg.png'), '-r300');
    % Save as PNG
    print(hFig1b, '-dpng', fullfile('..\figures\', 'mean_conn_mat_lap_mid.png'), '-r300');
    % Save as PNG
    print(hFig1c, '-dpng', fullfile('..\figures\', 'mean_conn_mat_lap_end.png'), '-r300');
end

%% TO DO!!!!
% network measures and AUC comparison
% check C:\Edgar\Dropbox\CIACYT\Projects\2022\fc\code\network_measures_script_2024





%% Perform statistical tests
H_lambda_HbO = nan([1, numel(threshold)]); P_lambda_HbO = nan([1, numel(threshold)]);
H_Ctot_pos_HbO = nan([1, numel(threshold)]); P_clust_coeff_HbO = nan([1, numel(threshold)]);
H_SWI_HbO = nan([1, numel(threshold)]); P_SWI_HbO = nan([1, numel(threshold)]);
nPerms = 1000;
for idxThresh = 1:numel(threshold)
    % Use permutation tests
    [P_lambda_HbO(idxThresh)] = permutationTest(lambda_HbO(design(:,2),idxThresh), lambda_HbO(design(:,1),idxThresh), nPerms);
    [P_clust_coeff_HbO(idxThresh)] = permutationTest(clust_coeff_HbO(design(:,2),idxThresh), clust_coeff_HbO(design(:,1),idxThresh), nPerms);
    [P_SWI_HbO(idxThresh)] = permutationTest(SWI_HbO(design(:,2),idxThresh), SWI_HbO(design(:,1),idxThresh), nPerms);
end
P_lambda_HbO = (P_lambda_HbO);
P_clust_coeff_HbO = (P_clust_coeff_HbO);
P_SWI_HbO = (P_SWI_HbO);

%% Plot P-values
alphaVal = 0.05;
figure(hFig2a);
% ------------------- Characteristic Path Length ------------------------
subplot(311); hold on
% yValue = nanmean(lambda_HbO) + nanstd(lambda_HbO);
yValue = 3.7*ones(size(threshold));
plot(threshold(P_lambda_HbO<=alphaVal), yValue (P_lambda_HbO<=alphaVal), 'k*');
axis tight
xlim([threshold(1) threshold(end)])
set(gca, 'FontSize', 16)
% legend({'Ctl', 'PD', 'p-Val<\alpha'}, 'Location', 'northwest');
% ------------------- Clustering Coefficient -------------------------------
subplot(312); hold on
% yValue = nanmean(lambda_HbO) + nanstd(lambda_HbO);
yValue = 1.1*ones(size(threshold));
plot(threshold((P_clust_coeff_HbO)<=alphaVal), yValue((P_clust_coeff_HbO)<=alphaVal), 'k*');
% axis tight
xlim([threshold(1) threshold(end)])
ylim([0 1.25])
legend({'Control', 'PD'}, 'Location', 'southeast');
set(gca, 'FontSize', 16)
% -------------- Small-world index ---------------------
subplot(313); hold on
% yValue = nanmean(lambda_HbO) + nanstd(lambda_HbO);
yValue = 3.5*ones(size(threshold));
plot(threshold((P_SWI_HbO)<=alphaVal), yValue ((P_SWI_HbO)<=alphaVal), 'k*');
% axis tight
xlim([threshold(1) threshold(end)])
ylim([0 3.8])
set(gca, 'FontSize', 16)
% legend({'Ctl', 'PD', 'p-Val<\alpha'}, 'Location', 'southeast');

% Specify window units
set(hFig2a, 'units', 'inches')
% Change figure and paper size
set(hFig2a, 'Position', [0.1 0.1 12 12])
set(hFig2a, 'PaperPosition', [0.1 0.1 12 12])

if printFigures
    % Save as PNG
    print(hFig2a, '-dpng', fullfile('..\figures\', 'small_world.png'), '-r300');
end

%% Compute AUC and statistical tests
AUC_lambda_HbO      = nan([nSubjects 1]);
AUC_clust_coeff_HbO = nan([nSubjects 1]);
AUC_SWI_HbO         = nan([nSubjects 1]);
for iSubjects = 1:nSubjects
    AUC_lambda_HbO(iSubjects) = trapz(threshold, lambda_HbO(iSubjects,:));
    AUC_clust_coeff_HbO(iSubjects) = trapz(threshold, clust_coeff_HbO(iSubjects,:));
    AUC_SWI_HbO(iSubjects) = trapz(threshold, SWI_HbO(iSubjects,:));
end
% Use permutation tests
P_AUC_lambda_HbO = permutationTest(AUC_lambda_HbO(design(:,2)), AUC_lambda_HbO(design(:,1)), nPerms);
P_AUC_clust_coeff_HbO = permutationTest(AUC_clust_coeff_HbO(design(:,2)), AUC_clust_coeff_HbO(design(:,1)), nPerms);
P_AUC_SWI_HbO = permutationTest(AUC_SWI_HbO(design(:,2)), AUC_SWI_HbO(design(:,1)), nPerms);

%% Plot AUC
markerStyle = ['o', 'x'];
markerColor = [0 0 1; 1 0 0]; % RGB for red, green, blue
markerSize = 50;
lineWidth = 2;
hFig3a = figure; set(hFig3a, 'color', 'w', 'Name', 'AUC char path length')
boxscatter([AUC_lambda_HbO(design(:,1)) AUC_lambda_HbO(design(:,2))], markerStyle, markerColor, markerSize, lineWidth);
% yLims = [0.05 0.91]; ylim(yLims)
title(sprintf('p=%0.4f', P_AUC_lambda_HbO))
ylabel('Characteristi Path Length')
hFig3b = figure; set(hFig3b, 'color', 'w', 'Name', 'Clustering Coefficient')
boxscatter([AUC_clust_coeff_HbO(design(:,2)) AUC_clust_coeff_HbO(design(:,1))], markerStyle, markerColor, markerSize, lineWidth);
% ylim(yLims)
title(sprintf('p=%0.4f', P_AUC_clust_coeff_HbO))
ylabel('Clustering Coefficient')
hFig3c = figure; set(hFig3c, 'color', 'w', 'Name', 'Small-World Index')
boxscatter([AUC_SWI_HbO(design(:,2)) AUC_SWI_HbO(design(:,1))], markerStyle, markerColor, markerSize, lineWidth);
% ylim(yLims)
title(sprintf('p=%0.4f', P_AUC_SWI_HbO))
ylabel('Small-World Index')

% Specify window units
set(hFig3a, 'units', 'inches')
% Change figure and paper size
set(hFig3a, 'Position', [0.1 0.1 4.5 4.5])
set(hFig3a, 'PaperPosition', [0.1 0.1 4.5 4.5])
% Specify window units
set(hFig3b, 'units', 'inches')
% Change figure and paper size
set(hFig3b, 'Position', [0.1 0.1 4.5 4.5])
set(hFig3b, 'PaperPosition', [0.1 0.1 4.5 4.5])
% Specify window units
set(hFig3c, 'units', 'inches')
% Change figure and paper size
set(hFig3c, 'Position', [0.1 0.1 4.5 4.5])
set(hFig3c, 'PaperPosition', [0.1 0.1 4.5 4.5])

if printFigures
    % Save as PNG
    print(hFig3a, '-dpng', fullfile('..\figures\', 'AUC_char_path_length.png'), '-r300');
    print(hFig3b, '-dpng', fullfile('..\figures\', 'AUC_clust_coeff.png'), '-r300');
    print(hFig3c, '-dpng', fullfile('..\figures\', 'AUC_swi.png'), '-r300');
end

%% NBS
d = 0.6;        % Cohen's d
t = sqrt(nSubjects) * d;

%% Compute correlations between lambda, clust_coeff, swi and mean connectivity
Mat_HbO_Ctrl = Mat_HbO(:,:,design(:,2));
Mat_HbO_Ctrl = reshape(Mat_HbO_Ctrl, [400 20]);
Mat_HbO_PD = Mat_HbO(:,:,design(:,1));
Mat_HbO_PD = reshape(Mat_HbO_PD, [400 20]);
mean_Conn_Ctrl = nanmean(Mat_HbO_Ctrl)';
mean_Conn_PD = nanmean(Mat_HbO_PD)';
[corr_lambda_Ctrl, corr_lambda_Ctrl_P] = corr(AUC_lambda_HbO(design(:,1)), mean_Conn_Ctrl);
[corr_lambda_PD, corr_lambda_PD_P] = corr(AUC_lambda_HbO(design(:,2)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_lambda_HbO(design(:,1)),1);
lambda_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_lambda_HbO(design(:,2)),1);
lambda_PD = polyval(p,mean_Conn_PD);
[corr_clust_Ctrl, corr_clust_Ctrl_P] = corr(AUC_clust_coeff_HbO(design(:,2)), mean_Conn_Ctrl);
[corr_clust_PD, corr_clust_PD_P] = corr(AUC_clust_coeff_HbO(design(:,1)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_clust_coeff_HbO(design(:,2)),1);
clust_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_clust_coeff_HbO(design(:,1)),1);
clust_PD = polyval(p,mean_Conn_PD);
[corr_swi_Ctrl, corr_swi_Ctrl_P] = corr(AUC_SWI_HbO(design(:,2)), mean_Conn_Ctrl);
[corr_swi_PD, corr_swi_PD_P] = corr(AUC_SWI_HbO(design(:,1)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_SWI_HbO(design(:,2)),1);
swi_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_SWI_HbO(design(:,1)),1);
swi_PD = polyval(p,mean_Conn_PD);


%% Plot correlations
% ----------------- Characteristic Path Length -----------------
hFig4a = figure; set(hFig4a, 'color', 'w', 'Name', 'Corr. Characteristic Path Length')
hold on
scatter(mean_Conn_Ctrl, AUC_lambda_HbO(design(:,1)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_lambda_HbO(design(:,2)), 'SizeData', markerSize, ...
            'Cdata', markerColor(2, :), 'Marker', markerStyle(2), ...
            'LineWidth', lineWidth);
title({sprintf('Control: r=%0.2f, p=%0.2f, mean FC=%0.2f', corr_lambda_Ctrl, corr_lambda_Ctrl_P, nanmean(mean_Conn_Ctrl)),...
    sprintf('PD: r=%0.2f, p=%0.2f, mean FC=%0.2f', corr_lambda_PD, corr_lambda_PD_P, nanmean(mean_Conn_PD))})
plot(mean_Conn_Ctrl, lambda_Ctrl, 'LineWidth', lineWidth, 'Color',markerColor(1, :))
plot(mean_Conn_PD, lambda_PD, 'LineWidth', lineWidth, 'Color',markerColor(2, :))
xlabel('Mean connectivity (z)')
ylabel('Characteristic Path Length')
set(gca, 'FontSize', 14)
% Specify window units
set(hFig4a, 'units', 'inches')
% Change figure and paper size
set(hFig4a, 'Position', [0.1 0.1 4.5 4.5])
set(hFig4a, 'PaperPosition', [0.1 0.1 4.5 4.5])
% ----------------- Clustering Coefficient -----------------
hFig4b = figure; set(hFig4b, 'color', 'w', 'Name', 'Corr. Clustering coefficient')
hold on
scatter(mean_Conn_Ctrl, AUC_clust_coeff_HbO(design(:,2)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_clust_coeff_HbO(design(:,1)), 'SizeData', markerSize, ...
            'Cdata', markerColor(2, :), 'Marker', markerStyle(2), ...
            'LineWidth', lineWidth);
title({sprintf('Control: r=%0.2f, *p=%0.2e, mean FC=%0.2f', corr_clust_Ctrl, corr_clust_Ctrl_P, nanmean(mean_Conn_Ctrl)),...
    sprintf('PD: r=%0.2f, *p=%0.2e, mean FC=%0.2f', corr_clust_PD, corr_clust_PD_P, nanmean(mean_Conn_PD))})
plot(mean_Conn_Ctrl, clust_Ctrl, 'LineWidth', lineWidth, 'Color',markerColor(1, :))
plot(mean_Conn_PD, clust_PD, 'LineWidth', lineWidth, 'Color',markerColor(2, :))
xlabel('Mean connectivity (z)')
ylabel('Clustering coefficient')
set(gca, 'FontSize', 14)
% Specify window units
set(hFig4b, 'units', 'inches')
% Change figure and paper size
set(hFig4b, 'Position', [0.1 0.1 4.5 4.5])
set(hFig4b, 'PaperPosition', [0.1 0.1 4.5 4.5])
% ----------------- Small-world index ----------------
hFig4c = figure; set(hFig4c, 'color', 'w', 'Name', 'Corr. Small World index')
hold on
scatter(mean_Conn_Ctrl, AUC_SWI_HbO(design(:,2)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_SWI_HbO(design(:,1)), 'SizeData', markerSize, ...
            'Cdata', markerColor(2, :), 'Marker', markerStyle(2), ...
            'LineWidth', lineWidth);
title({sprintf('Control: r=%0.2f, *p=%0.2e, mean FC=%0.2f', corr_swi_Ctrl, corr_swi_Ctrl_P, nanmean(mean_Conn_Ctrl)),...
    sprintf('PD: r=%0.2f, *p=%0.2e, mean FC=%0.2f', corr_swi_PD, corr_swi_PD_P, nanmean(mean_Conn_PD))})
plot(mean_Conn_Ctrl, swi_Ctrl, 'LineWidth', lineWidth, 'Color',markerColor(1, :))
plot(mean_Conn_PD, swi_PD, 'LineWidth', lineWidth, 'Color',markerColor(2, :))
xlabel('Mean connectivity (z)')
ylabel('Clustering coefficient')
set(gca, 'FontSize', 14)
% Specify window units
set(hFig4c, 'units', 'inches')
% Change figure and paper size
set(hFig4c, 'Position', [0.1 0.1 4.5 4.5])
set(hFig4c, 'PaperPosition', [0.1 0.1 4.5 4.5])

%% Print correlations
if printFigures
    % Save as PNG
    print(hFig4a, '-dpng', fullfile('..\figures\', 'corr_char_path_length.png'), '-r300');
    print(hFig4b, '-dpng', fullfile('..\figures\', 'corr_clust_coeff.png'), '-r300');
    print(hFig4c, '-dpng', fullfile('..\figures\', 'corr_swi.png'), '-r300');
end

% EOF
