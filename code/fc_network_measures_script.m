%% Load data
clear; close all; clc;
% Load data

%% Prepare data
% Load data
load('..\data\lap_all_significant_connections_HbO.mat');
rng(42); % For repeatability
printFigures = false;
symmColorZero = true;

%% Prepare data
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];
% meanMat22 is your 22×22 mean connectivity matrix.
% meanMatBegLap = meanMatFromCell(beginLap.zMatFDR, channelIdx);
% meanMatMidLap = meanMatFromCell(midLap.zMatFDR, channelIdx);
% meanMatEndLap = meanMatFromCell(endLap.zMatFDR, channelIdx);
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
designMat = false([size(beginLap.zMatFDR,1)+size(midLap.zMatFDR,1)+size(endLap.zMatFDR,1), 3]);
designMat(1:30, 1) = true;
designMat(31:60, 2) = true;
designMat(61:end, 3) = true;

%% tensor initialization
% Initialize tensor to hold the 22x22x90 matrices
Mat_HbO = zeros(22, 22, 90);

% Combine all 3 cells into a single 90-element cell array
allMats = [beginLap.zMatFDR; midLap.zMatFDR; endLap.zMatFDR]; % 90x1 cell

% Loop through each matrix, extract the desired channels, and fill the tensor
for k = 1:90
    mat24 = allMats{k};                % 24x24 matrix
    mat22 = mat24(channelIdx, channelIdx); % 22x22 submatrix
    Mat_HbO(:, :, k) = mat22;
end


%% Compute Characteristic path length (lambda), clustering coefficient
% Sparsity
threshold = 0.1:0.01:0.34;
% number of scans
nSubjects = size(Mat_HbO, 3);
% Initialize variables
lambda_HbO = nan([nSubjects, numel(threshold)]);
efficiency_HbO = nan([nSubjects, numel(threshold)]);
clust_coeff_HbO = nan([nSubjects, numel(threshold)]);
SWI_HbO = nan([nSubjects, numel(threshold)]);
for idxThresh = 1:numel(threshold)
    for iSubjects=1:nSubjects
        W = squeeze(Mat_HbO(:,:,iSubjects));   % resting state (HbO)
        % -------------- Characteristic Path Length ---------------------
        Wthresh = threshold_proportional(W, threshold(idxThresh));  % Applying proportional threshold
        L = weight_conversion(Wthresh, 'lengths');            % taking the inverse of matrix elements
        D = distance_wei(L);                            % applying Dijkstra's algorithm
        % compute characteristic path length and global efficiency
        [lambda_HbO(iSubjects, idxThresh), efficiency_HbO(iSubjects, idxThresh)] = charpath(D, 0, 0);
        % Compute Modularity
        % [~, Q_L_hbo(iSubjects, idx)] = modularity_und(W);
        % -------------- Clustering coefficient ---------------------
        % Compute clustering coefficient
        clust_coeff_HbO(iSubjects, idxThresh) = nanmean(clustering_coef_bu(Wthresh));
        % -------------- Small-world index ---------------------
        SWI_HbO(iSubjects, idxThresh) = swi(Wthresh);
    end
end

%% Plot Characteristic path length and global efficiency as a function of threshold
controlColor = [0 0 1];
hFig2a = figure; set(hFig2a, 'color', 'w', 'Name', 'smallWorld')
% ------------------- Characteristic Path Length ------------------------
subplot(311); hold on
% small-world properties (mean +- s.d.) 
hCtl = shadedErrorBar(threshold, lambda_HbO(designMat(:,1),:),{@nanmean,@nanstd},'-b',1);
% hCtl = shadedErrorBarEG(timeVector, 1e6*squeeze(HbOchannels.R(~isParkinson,iChannels,:)),{@mean, stderror},'-', 1);
% hCtl.mainLine.Color = controlColor;
% hCtl.edge(1).Color = controlColor;
% hCtl.edge(2).Color = controlColor;
% hCtl.patch.FaceColor = controlColor;
hCtl.mainLine.LineWidth = 2;
shadedErrorBar(threshold,lambda_HbO(designMat(:,2),:),{@nanmean,@nanstd},'-r',1);
plot(threshold, nanmean(lambda_HbO(designMat(:,1),:)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(lambda_HbO(designMat(:,2),:)), 'r-', 'LineWidth', 2);
% title('\lambda')
xlabel('Sparsity'); ylabel('Characteristic Path Length')
% ------------------- Clustering Coefficient -------------------------------
subplot(312); hold on
% small-world properties (mean +- s.d.) 
shadedErrorBar(threshold, clust_coeff_HbO(designMat(:,2),:),{@nanmean,@nanstd},'-b',1);
shadedErrorBar(threshold,clust_coeff_HbO(designMat(:,1),:),{@nanmean,@nanstd},'-r',1);
plot(threshold, nanmean(clust_coeff_HbO(designMat(:,2),:)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(clust_coeff_HbO(designMat(:,1),:)), 'r-', 'LineWidth', 2);
% title('CC')
xlabel('Sparsity'); ylabel('Clustering Coefficient')
% -------------- Small-world index ---------------------
subplot(313); hold on
% small-world properties (mean +- s.d.) 
shadedErrorBar(threshold, SWI_HbO(designMat(:,2),:),{@nanmean,@nanstd},'-b',1);
shadedErrorBar(threshold,SWI_HbO(designMat(:,1),:),{@nanmean,@nanstd},'-r',1);
plot(threshold, nanmean(SWI_HbO(designMat(:,2),:)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(SWI_HbO(designMat(:,1),:)), 'r-', 'LineWidth', 2);
% title('SWI')
xlabel('Sparsity'); ylabel('Small-world index')


%% Perform statistical tests
H_lambda_HbO = nan([1, numel(threshold)]); P_lambda_HbO = nan([1, numel(threshold)]);
H_Ctot_pos_HbO = nan([1, numel(threshold)]); P_clust_coeff_HbO = nan([1, numel(threshold)]);
H_SWI_HbO = nan([1, numel(threshold)]); P_SWI_HbO = nan([1, numel(threshold)]);
nPerms = 1000;
for idxThresh = 1:numel(threshold)
    % Use permutation tests
    [P_lambda_HbO(idxThresh)] = permutationTest(lambda_HbO(designMat(:,2),idxThresh), lambda_HbO(designMat(:,1),idxThresh), nPerms);
    [P_clust_coeff_HbO(idxThresh)] = permutationTest(clust_coeff_HbO(designMat(:,2),idxThresh), clust_coeff_HbO(designMat(:,1),idxThresh), nPerms);
    [P_SWI_HbO(idxThresh)] = permutationTest(SWI_HbO(designMat(:,2),idxThresh), SWI_HbO(designMat(:,1),idxThresh), nPerms);
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
P_AUC_lambda_HbO = permutationTest(AUC_lambda_HbO(designMat(:,2)), AUC_lambda_HbO(designMat(:,1)), nPerms);
P_AUC_clust_coeff_HbO = permutationTest(AUC_clust_coeff_HbO(designMat(:,2)), AUC_clust_coeff_HbO(designMat(:,1)), nPerms);
P_AUC_SWI_HbO = permutationTest(AUC_SWI_HbO(designMat(:,2)), AUC_SWI_HbO(designMat(:,1)), nPerms);

%% Plot AUC
markerStyle = ['o', 'x'];
markerColor = [0 0 1; 1 0 0]; % RGB for red, green, blue
markerSize = 50;
lineWidth = 2;
hFig3a = figure; set(hFig3a, 'color', 'w', 'Name', 'AUC char path length')
boxscatter([AUC_lambda_HbO(designMat(:,1)) AUC_lambda_HbO(designMat(:,2))], markerStyle, markerColor, markerSize, lineWidth);
% yLims = [0.05 0.91]; ylim(yLims)
title(sprintf('p=%0.4f', P_AUC_lambda_HbO))
ylabel('Characteristi Path Length')
hFig3b = figure; set(hFig3b, 'color', 'w', 'Name', 'Clustering Coefficient')
boxscatter([AUC_clust_coeff_HbO(designMat(:,2)) AUC_clust_coeff_HbO(designMat(:,1))], markerStyle, markerColor, markerSize, lineWidth);
% ylim(yLims)
title(sprintf('p=%0.4f', P_AUC_clust_coeff_HbO))
ylabel('Clustering Coefficient')
hFig3c = figure; set(hFig3c, 'color', 'w', 'Name', 'Small-World Index')
boxscatter([AUC_SWI_HbO(designMat(:,2)) AUC_SWI_HbO(designMat(:,1))], markerStyle, markerColor, markerSize, lineWidth);
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
Mat_HbO_Ctrl = Mat_HbO(:,:,designMat(:,2));
Mat_HbO_Ctrl = reshape(Mat_HbO_Ctrl, [400 20]);
Mat_HbO_PD = Mat_HbO(:,:,designMat(:,1));
Mat_HbO_PD = reshape(Mat_HbO_PD, [400 20]);
mean_Conn_Ctrl = nanmean(Mat_HbO_Ctrl)';
mean_Conn_PD = nanmean(Mat_HbO_PD)';
[corr_lambda_Ctrl, corr_lambda_Ctrl_P] = corr(AUC_lambda_HbO(designMat(:,1)), mean_Conn_Ctrl);
[corr_lambda_PD, corr_lambda_PD_P] = corr(AUC_lambda_HbO(designMat(:,2)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_lambda_HbO(designMat(:,1)),1);
lambda_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_lambda_HbO(designMat(:,2)),1);
lambda_PD = polyval(p,mean_Conn_PD);
[corr_clust_Ctrl, corr_clust_Ctrl_P] = corr(AUC_clust_coeff_HbO(designMat(:,2)), mean_Conn_Ctrl);
[corr_clust_PD, corr_clust_PD_P] = corr(AUC_clust_coeff_HbO(designMat(:,1)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_clust_coeff_HbO(designMat(:,2)),1);
clust_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_clust_coeff_HbO(designMat(:,1)),1);
clust_PD = polyval(p,mean_Conn_PD);
[corr_swi_Ctrl, corr_swi_Ctrl_P] = corr(AUC_SWI_HbO(designMat(:,2)), mean_Conn_Ctrl);
[corr_swi_PD, corr_swi_PD_P] = corr(AUC_SWI_HbO(designMat(:,1)), mean_Conn_PD);
p = polyfit(mean_Conn_Ctrl, AUC_SWI_HbO(designMat(:,2)),1);
swi_Ctrl = polyval(p,mean_Conn_Ctrl);
p = polyfit(mean_Conn_PD, AUC_SWI_HbO(designMat(:,1)),1);
swi_PD = polyval(p,mean_Conn_PD);


%% Plot correlations
% ----------------- Characteristic Path Length -----------------
hFig4a = figure; set(hFig4a, 'color', 'w', 'Name', 'Corr. Characteristic Path Length')
hold on
scatter(mean_Conn_Ctrl, AUC_lambda_HbO(designMat(:,1)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_lambda_HbO(designMat(:,2)), 'SizeData', markerSize, ...
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
scatter(mean_Conn_Ctrl, AUC_clust_coeff_HbO(designMat(:,2)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_clust_coeff_HbO(designMat(:,1)), 'SizeData', markerSize, ...
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
scatter(mean_Conn_Ctrl, AUC_SWI_HbO(designMat(:,2)), 'SizeData', markerSize, ...
            'Cdata', markerColor(1, :), 'Marker', markerStyle(1), ...
            'LineWidth', lineWidth);
scatter(mean_Conn_PD, AUC_SWI_HbO(designMat(:,1)), 'SizeData', markerSize, ...
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
