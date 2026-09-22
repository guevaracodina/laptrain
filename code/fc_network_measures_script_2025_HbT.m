%% Load data
clear; close all; clc;
% Load data

%% Prepare data
% Load data
load('..\data\lap_all_significant_connections_HbT.mat');
rng(42); % For repeatability
printFigures = false;
symmColorZero = true;

%% Prepare data
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];
% meanMat22 is your 22×22 mean connectivity matrix (after selecting channels).
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
% Design matrix for 3 groups: Start, Middle, End (each group has 30 measures)
designMat = false(size(beginLap.zMatFDR,1) + size(midLap.zMatFDR,1) + size(endLap.zMatFDR,1), 3);
designMat(1:30, 1)   = true;  % Start group indices (1–30)
designMat(31:60, 2)  = true;  % Middle group indices (31–60)
designMat(61:end, 3) = true;  % End group indices (61–90)

%% Tensor initialization for connectivity matrices
% Initialize tensor to hold the 22x22x90 matrices (90 subjects total)
Mat_HbT = zeros(22, 22, sum(designMat(:)));
% Combine all 3 group matrices into a single 90-element cell array
allMats = [beginLap.zMatFDR; midLap.zMatFDR; endLap.zMatFDR];  % 90x1 cell
% Loop through each matrix, extract desired channels, and fill the tensor
for k = 1:numel(allMats)
    mat24 = allMats{k};                     % original 24x24 matrix
    mat22 = mat24(channelIdx, channelIdx);  % 22x22 submatrix for selected channels
    Mat_HbT(:, :, k) = mat22;
end

%% Compute network metrics (Characteristic Path Length, Clustering Coefficient, Small-World Index)
threshold = 0.1:0.01:0.34;
nSubjects = size(Mat_HbT, 3);
% Preallocate matrices: [subjects x thresholds]
lambda_HbT      = nan(nSubjects, numel(threshold));
efficiency_HbT  = nan(nSubjects, numel(threshold));
clust_coeff_HbT = nan(nSubjects, numel(threshold));
SWI_HbT         = nan(nSubjects, numel(threshold));
for idxThresh = 1:numel(threshold)
    for iSub = 1:nSubjects
        W = Mat_HbT(:, :, iSub);                           % connectivity matrix for subject
        Wthresh = threshold_proportional(W, threshold(idxThresh));  % apply proportional threshold
        L = weight_conversion(Wthresh, 'lengths');         % invert weights to lengths
        D = distance_wei(L);                               % compute distance matrix for lengths
        % Characteristic path length (lambda) and global efficiency
        [lambda_HbT(iSub, idxThresh), efficiency_HbT(iSub, idxThresh)] = charpath(D, 0, 0);
        % Clustering coefficient (average)
        clust_coeff_HbT(iSub, idxThresh) = nanmean(clustering_coef_bu(Wthresh));
        % Small-world index
        SWI_HbT(iSub, idxThresh) = swi(Wthresh);
    end
end

%% Plot metrics vs. threshold for 3 groups
% Define colors for groups: Start=blue, Middle=green, End=red
startColor = [0 0 1]; midColor = [0 1 0]; endColor = [1 0 0];
hFig2a = figure; set(hFig2a, 'color', 'w', 'Name', 'smallWorld');
% --- Characteristic Path Length ---
subplot(311); hold on;
shadedErrorBar(threshold, lambda_HbT(designMat(:,1), :), {@nanmean, @nanstd}, '-b', 1);
shadedErrorBar(threshold, lambda_HbT(designMat(:,2), :), {@nanmean, @nanstd}, '-g', 1);
shadedErrorBar(threshold, lambda_HbT(designMat(:,3), :), {@nanmean, @nanstd}, '-r', 1);
% Plot mean lines for each group (thicker lines for clarity)
plot(threshold, nanmean(lambda_HbT(designMat(:,1), :)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(lambda_HbT(designMat(:,2), :)), 'g-', 'LineWidth', 2);
plot(threshold, nanmean(lambda_HbT(designMat(:,3), :)), 'r-', 'LineWidth', 2);
xlabel('Sparsity'); ylabel('Characteristic Path Length');
% --- Clustering Coefficient ---
subplot(312); hold on;
shadedErrorBar(threshold, clust_coeff_HbT(designMat(:,1), :), {@nanmean, @nanstd}, '-b', 1);
shadedErrorBar(threshold, clust_coeff_HbT(designMat(:,2), :), {@nanmean, @nanstd}, '-g', 1);
shadedErrorBar(threshold, clust_coeff_HbT(designMat(:,3), :), {@nanmean, @nanstd}, '-r', 1);
plot(threshold, nanmean(clust_coeff_HbT(designMat(:,1), :)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(clust_coeff_HbT(designMat(:,2), :)), 'g-', 'LineWidth', 2);
plot(threshold, nanmean(clust_coeff_HbT(designMat(:,3), :)), 'r-', 'LineWidth', 2);
xlabel('Sparsity'); ylabel('Clustering Coefficient');
% --- Small-World Index ---
subplot(313); hold on;
shadedErrorBar(threshold, SWI_HbT(designMat(:,1), :), {@nanmean, @nanstd}, '-b', 1);
shadedErrorBar(threshold, SWI_HbT(designMat(:,2), :), {@nanmean, @nanstd}, '-g', 1);
shadedErrorBar(threshold, SWI_HbT(designMat(:,3), :), {@nanmean, @nanstd}, '-r', 1);
plot(threshold, nanmean(SWI_HbT(designMat(:,1), :)), 'b-', 'LineWidth', 2);
plot(threshold, nanmean(SWI_HbT(designMat(:,2), :)), 'g-', 'LineWidth', 2);
plot(threshold, nanmean(SWI_HbT(designMat(:,3), :)), 'r-', 'LineWidth', 2);
xlabel('Sparsity'); ylabel('Small-World Index');
% Add legend (using one subplot for clarity)
subplot(312);
legend({'Start', 'Middle', 'End'}, 'Location', 'best');

%% Statistical tests at each threshold (3-group and pairwise post-hoc)
alphaVal = 0.05;
nPerms   = 1000;
P_lambda_HbT      = nan(1, numel(threshold));
P_clust_coeff_HbT = nan(1, numel(threshold));
P_SWI_HbT         = nan(1, numel(threshold));
% Initialize matrices for pairwise p-values (post-hoc) at each threshold
P_lambda_HbT_12 = nan(1, numel(threshold));  % Start vs Middle
P_lambda_HbT_13 = nan(1, numel(threshold));  % Start vs End
P_lambda_HbT_23 = nan(1, numel(threshold));  % Middle vs End
P_clust_coeff_HbT_12 = nan(1, numel(threshold));
P_clust_coeff_HbT_13 = nan(1, numel(threshold));
P_clust_coeff_HbT_23 = nan(1, numel(threshold));
P_SWI_HbT_12 = nan(1, numel(threshold));
P_SWI_HbT_13 = nan(1, numel(threshold));
P_SWI_HbT_23 = nan(1, numel(threshold));
for idxThresh = 1:numel(threshold)
    % Permutation test across 3 groups (global null hypothesis)
    P_lambda_HbT(idxThresh) = permutationTest3( ...
        lambda_HbT(designMat(:,1), idxThresh), ...
        lambda_HbT(designMat(:,2), idxThresh), ...
        lambda_HbT(designMat(:,3), idxThresh), nPerms);
    P_clust_coeff_HbT(idxThresh) = permutationTest3( ...
        clust_coeff_HbT(designMat(:,1), idxThresh), ...
        clust_coeff_HbT(designMat(:,2), idxThresh), ...
        clust_coeff_HbT(designMat(:,3), idxThresh), nPerms);
    P_SWI_HbT(idxThresh) = permutationTest3( ...
        SWI_HbT(designMat(:,1), idxThresh), ...
        SWI_HbT(designMat(:,2), idxThresh), ...
        SWI_HbT(designMat(:,3), idxThresh), nPerms);
    % If global test is significant, perform pairwise permutation tests (post-hoc)
    if P_lambda_HbT(idxThresh) <= alphaVal
        P_lambda_HbT_12(idxThresh) = permutationTest(lambda_HbT(designMat(:,1), idxThresh), lambda_HbT(designMat(:,2), idxThresh), nPerms);
        P_lambda_HbT_13(idxThresh) = permutationTest(lambda_HbT(designMat(:,1), idxThresh), lambda_HbT(designMat(:,3), idxThresh), nPerms);
        P_lambda_HbT_23(idxThresh) = permutationTest(lambda_HbT(designMat(:,2), idxThresh), lambda_HbT(designMat(:,3), idxThresh), nPerms);
    end
    if P_clust_coeff_HbT(idxThresh) <= alphaVal
        P_clust_coeff_HbT_12(idxThresh) = permutationTest(clust_coeff_HbT(designMat(:,1), idxThresh), clust_coeff_HbT(designMat(:,2), idxThresh), nPerms);
        P_clust_coeff_HbT_13(idxThresh) = permutationTest(clust_coeff_HbT(designMat(:,1), idxThresh), clust_coeff_HbT(designMat(:,3), idxThresh), nPerms);
        P_clust_coeff_HbT_23(idxThresh) = permutationTest(clust_coeff_HbT(designMat(:,2), idxThresh), clust_coeff_HbT(designMat(:,3), idxThresh), nPerms);
    end
    if P_SWI_HbT(idxThresh) <= alphaVal
        P_SWI_HbT_12(idxThresh) = permutationTest(SWI_HbT(designMat(:,1), idxThresh), SWI_HbT(designMat(:,2), idxThresh), nPerms);
        P_SWI_HbT_13(idxThresh) = permutationTest(SWI_HbT(designMat(:,1), idxThresh), SWI_HbT(designMat(:,3), idxThresh), nPerms);
        P_SWI_HbT_23(idxThresh) = permutationTest(SWI_HbT(designMat(:,2), idxThresh), SWI_HbT(designMat(:,3), idxThresh), nPerms);
    end
end

%% Mark significant thresholds on plots (overall 3-group test)
figure(hFig2a);
% Characteristic Path Length subplot
subplot(311); hold on;
yLine = 3.7 * ones(size(threshold));  % baseline for star markers
plot(threshold(P_lambda_HbT <= alphaVal), yLine(P_lambda_HbT <= alphaVal), 'k*');
axis tight;  % tighten y-axis to data
xlim([threshold(1) threshold(end)]);
set(gca, 'FontSize', 16);
% Clustering Coefficient subplot
subplot(312); hold on;
yLine = 1.1 * ones(size(threshold));
plot(threshold(P_clust_coeff_HbT <= alphaVal), yLine(P_clust_coeff_HbT <= alphaVal), 'k*');
xlim([threshold(1) threshold(end)]);
ylim([0 1.25]);  % include baseline at 0
set(gca, 'FontSize', 16);
% Small-World Index subplot
subplot(313); hold on;
yLine = 3.5 * ones(size(threshold));
plot(threshold(P_SWI_HbT <= alphaVal), yLine(P_SWI_HbT <= alphaVal), 'k*');
xlim([threshold(1) threshold(end)]);
ylim([0 3.8]);
set(gca, 'FontSize', 16);
% (Stars * indicate thresholds where there is a significant overall difference among the 3 groups at p < 0.05)

%% Compute AUC for each metric per subject and statistical tests on AUC
AUC_lambda_HbT      = nan(nSubjects, 1);
AUC_clust_coeff_HbT = nan(nSubjects, 1);
AUC_SWI_HbT         = nan(nSubjects, 1);
for iSub = 1:nSubjects
    AUC_lambda_HbT(iSub)      = trapz(threshold, lambda_HbT(iSub, :));
    AUC_clust_coeff_HbT(iSub) = trapz(threshold, clust_coeff_HbT(iSub, :));
    AUC_SWI_HbT(iSub)         = trapz(threshold, SWI_HbT(iSub, :));
end
% Permutation test across 3 groups for AUC (global p-value)
P_AUC_lambda_HbT      = permutationTest3(AUC_lambda_HbT(designMat(:,1)), AUC_lambda_HbT(designMat(:,2)), AUC_lambda_HbT(designMat(:,3)), nPerms);
P_AUC_clust_coeff_HbT = permutationTest3(AUC_clust_coeff_HbT(designMat(:,1)), AUC_clust_coeff_HbT(designMat(:,2)), AUC_clust_coeff_HbT(designMat(:,3)), nPerms);
P_AUC_SWI_HbT         = permutationTest3(AUC_SWI_HbT(designMat(:,1)), AUC_SWI_HbT(designMat(:,2)), AUC_SWI_HbT(designMat(:,3)), nPerms);
% Pairwise permutation tests for AUC (post-hoc comparisons)
P_AUC_lambda_12 = permutationTest(AUC_lambda_HbT(designMat(:,1)), AUC_lambda_HbT(designMat(:,2)), nPerms);
P_AUC_lambda_13 = permutationTest(AUC_lambda_HbT(designMat(:,1)), AUC_lambda_HbT(designMat(:,3)), nPerms);
P_AUC_lambda_23 = permutationTest(AUC_lambda_HbT(designMat(:,2)), AUC_lambda_HbT(designMat(:,3)), nPerms);
P_AUC_clust_12  = permutationTest(AUC_clust_coeff_HbT(designMat(:,1)), AUC_clust_coeff_HbT(designMat(:,2)), nPerms);
P_AUC_clust_13  = permutationTest(AUC_clust_coeff_HbT(designMat(:,1)), AUC_clust_coeff_HbT(designMat(:,3)), nPerms);
P_AUC_clust_23  = permutationTest(AUC_clust_coeff_HbT(designMat(:,2)), AUC_clust_coeff_HbT(designMat(:,3)), nPerms);
P_AUC_SWI_12    = permutationTest(AUC_SWI_HbT(designMat(:,1)), AUC_SWI_HbT(designMat(:,2)), nPerms);
P_AUC_SWI_13    = permutationTest(AUC_SWI_HbT(designMat(:,1)), AUC_SWI_HbT(designMat(:,3)), nPerms);
P_AUC_SWI_23    = permutationTest(AUC_SWI_HbT(designMat(:,2)), AUC_SWI_HbT(designMat(:,3)), nPerms);

%% Plot AUC comparisons for each metric (3 groups)
markerStyle = ['o', 'x', '^'];
markerColor = [0 0 1; 0 1 0; 1 0 0];  % [Blue; Green; Red]
markerSize  = 50;
lineWidth   = 2;
% Characteristic Path Length AUC plot
hFig3a = figure; set(hFig3a, 'color', 'w', 'Name', 'AUC_CharPathLength');
boxscatter([AUC_lambda_HbT(designMat(:,1)), AUC_lambda_HbT(designMat(:,2)), AUC_lambda_HbT(designMat(:,3))], ...
           markerStyle, markerColor, markerSize, lineWidth);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Start', 'Middle', 'End'});
ylabel('Characteristic Path Length (AUC)');
title({sprintf('Global p = %0.4f', P_AUC_lambda_HbT), ...
       sprintf('Post-hoc p: Start vs Middle = %0.4f, Start vs End = %0.4f, Middle vs End = %0.4f', ...
               P_AUC_lambda_12, P_AUC_lambda_13, P_AUC_lambda_23)});
% Clustering Coefficient AUC plot
hFig3b = figure; set(hFig3b, 'color', 'w', 'Name', 'AUC_ClustCoeff');
boxscatter([AUC_clust_coeff_HbT(designMat(:,1)), AUC_clust_coeff_HbT(designMat(:,2)), AUC_clust_coeff_HbT(designMat(:,3))], ...
           markerStyle, markerColor, markerSize, lineWidth);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Start', 'Middle', 'End'});
ylabel('Clustering Coefficient (AUC)');
title({sprintf('Global p = %0.4f', P_AUC_clust_coeff_HbT), ...
       sprintf('Post-hoc p: Start vs Middle = %0.4f, Start vs End = %0.4f, Middle vs End = %0.4f', ...
               P_AUC_clust_12, P_AUC_clust_13, P_AUC_clust_23)});
% Small-World Index AUC plot
hFig3c = figure; set(hFig3c, 'color', 'w', 'Name', 'AUC_SWI');
boxscatter([AUC_SWI_HbT(designMat(:,1)), AUC_SWI_HbT(designMat(:,2)), AUC_SWI_HbT(designMat(:,3))], ...
           markerStyle, markerColor, markerSize, lineWidth);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Start', 'Middle', 'End'});
ylabel('Small-World Index (AUC)');
title({sprintf('Global p = %0.4f', P_AUC_SWI_HbT), ...
       sprintf('Post-hoc p: Start vs Middle = %0.4f, Start vs End = %0.4f, Middle vs End = %0.4f', ...
               P_AUC_SWI_12, P_AUC_SWI_13, P_AUC_SWI_23)});
% (Each figure's title shows the global p-value and pairwise p-values for that metric's AUC)

% Adjust figure sizes if needed (for printing)
set([hFig3a, hFig3b, hFig3c], 'Units', 'inches', 'Position', [0.1 0.1 5 5], 'PaperPosition', [0.1 0.1 5 5]);
if printFigures
    print(hFig3a, '-dpng', fullfile('..\figures', 'AUC_char_path_length.png'), '-r300');
    print(hFig3b, '-dpng', fullfile('..\figures', 'AUC_clust_coeff.png'), '-r300');
    print(hFig3c, '-dpng', fullfile('..\figures', 'AUC_swi.png'), '-r300');
end



%% Correlation analysis between AUC metrics and mean connectivity (for each group)
% Compute connectivity matrices for each group (considering significant connections)
Mat_HbT_Start = Mat_HbT(:, :, designMat(:,1));  % 22x22x30 for Start group
Mat_HbT_Mid   = Mat_HbT(:, :, designMat(:,2));  % 22x22x30 for Middle group
Mat_HbT_End   = Mat_HbT(:, :, designMat(:,3));  % 22x22x30 for End group

% Reshape to 2D (edges x subjects) – assuming 484 edges (22x22 matrix) per subject
Mat_HbT_Start = reshape(Mat_HbT_Start, 484, sum(designMat(:,1)));  
Mat_HbT_Mid   = reshape(Mat_HbT_Mid,   484, sum(designMat(:,2)));
Mat_HbT_End   = reshape(Mat_HbT_End,   484, sum(designMat(:,3)));

% **Remove outlier values** in each group's matrix using isoutlier with 'mean' method (3 SD from mean)
% TF_start = isoutlier(Mat_HbT_Start, 'median');   % logical matrix of outliers for Start group (column-wise by subject):contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}
% Mat_HbT_Start(TF_start) = NaN;                 % set outlier entries to NaN to exclude them from mean calculation
% TF_mid   = isoutlier(Mat_HbT_Mid, 'median');     % outliers for Middle group
% Mat_HbT_Mid(TF_mid)     = NaN;
% TF_end   = isoutlier(Mat_HbT_End, 'median');     % outliers for End group
% Mat_HbT_End(TF_end)     = NaN;

% Compute mean connectivity for each subject in each group (ignoring NaNs from outliers)
mean_Conn_Start = nanmedian(Mat_HbT_Start)';  % 30x1 vector (average FC per Start subject, excluding outlier values)
mean_Conn_Mid   = nanmedian(Mat_HbT_Mid)';    % 30x1 vector (Middle group)
mean_Conn_End   = nanmedian(Mat_HbT_End)';    % 30x1 vector (End group)

% Removing outliers
TF_start = isoutlier(mean_Conn_Start, 'median');
mean_Conn_Start = mean_Conn_Start(~TF_start);
TF_mid   = isoutlier(mean_Conn_Mid, 'median');     % outliers for Middle group
mean_Conn_Mid = mean_Conn_Mid(~TF_mid);
TF_end   = isoutlier(mean_Conn_End, 'median');     % outliers for End group
mean_Conn_End = mean_Conn_End(~TF_end);
designMat(1:30,1) = ~TF_start;
designMat(31:60,2) = ~TF_mid;
designMat(61:90,3) = ~TF_end;


% Correlations between AUC metrics and mean connectivity within each group
[corr_lambda_Start, corr_lambda_Start_P] = corr(AUC_lambda_HbT(designMat(:,1)), mean_Conn_Start);
[corr_lambda_Mid,   corr_lambda_Mid_P]   = corr(AUC_lambda_HbT(designMat(:,2)), mean_Conn_Mid);
[corr_lambda_End,   corr_lambda_End_P]   = corr(AUC_lambda_HbT(designMat(:,3)), mean_Conn_End);
[corr_clust_Start,  corr_clust_Start_P]  = corr(AUC_clust_coeff_HbT(designMat(:,1)), mean_Conn_Start);
[corr_clust_Mid,    corr_clust_Mid_P]    = corr(AUC_clust_coeff_HbT(designMat(:,2)), mean_Conn_Mid);
[corr_clust_End,    corr_clust_End_P]    = corr(AUC_clust_coeff_HbT(designMat(:,3)), mean_Conn_End);
[corr_swi_Start,    corr_swi_Start_P]    = corr(AUC_SWI_HbT(designMat(:,1)), mean_Conn_Start);
[corr_swi_Mid,      corr_swi_Mid_P]      = corr(AUC_SWI_HbT(designMat(:,2)), mean_Conn_Mid);
[corr_swi_End,      corr_swi_End_P]      = corr(AUC_SWI_HbT(designMat(:,3)), mean_Conn_End);
% Prepare regression line fits for scatter plots
xRange_Start = linspace(min(mean_Conn_Start), max(mean_Conn_Start), 100);
xRange_Mid   = linspace(min(mean_Conn_Mid),   max(mean_Conn_Mid),   100);
xRange_End   = linspace(min(mean_Conn_End),   max(mean_Conn_End),   100);
lambda_Start_fit = polyval(polyfit(mean_Conn_Start, AUC_lambda_HbT(designMat(:,1)), 1), xRange_Start);
lambda_Mid_fit   = polyval(polyfit(mean_Conn_Mid,   AUC_lambda_HbT(designMat(:,2)), 1), xRange_Mid);
lambda_End_fit   = polyval(polyfit(mean_Conn_End,   AUC_lambda_HbT(designMat(:,3)), 1), xRange_End);
clust_Start_fit  = polyval(polyfit(mean_Conn_Start, AUC_clust_coeff_HbT(designMat(:,1)), 1), xRange_Start);
clust_Mid_fit    = polyval(polyfit(mean_Conn_Mid,   AUC_clust_coeff_HbT(designMat(:,2)), 1), xRange_Mid);
clust_End_fit    = polyval(polyfit(mean_Conn_End,   AUC_clust_coeff_HbT(designMat(:,3)), 1), xRange_End);
swi_Start_fit    = polyval(polyfit(mean_Conn_Start, AUC_SWI_HbT(designMat(:,1)), 1), xRange_Start);
swi_Mid_fit      = polyval(polyfit(mean_Conn_Mid,   AUC_SWI_HbT(designMat(:,2)), 1), xRange_Mid);
swi_End_fit      = polyval(polyfit(mean_Conn_End,   AUC_SWI_HbT(designMat(:,3)), 1), xRange_End);
% Plot correlation: Characteristic Path Length AUC vs mean connectivity
hFig4a = figure; set(hFig4a, 'color', 'w', 'Name', 'Corr_CharPathLength');
hold on;
scatter(mean_Conn_Start, AUC_lambda_HbT(designMat(:,1)), markerSize, 'MarkerEdgeColor', startColor, 'MarkerFaceColor', startColor, 'Marker', 'o', 'LineWidth', lineWidth);
scatter(mean_Conn_Mid,   AUC_lambda_HbT(designMat(:,2)), markerSize, 'MarkerEdgeColor', midColor, 'Marker', 'x', 'LineWidth', lineWidth);
scatter(mean_Conn_End,   AUC_lambda_HbT(designMat(:,3)), markerSize, 'MarkerEdgeColor', endColor, 'MarkerFaceColor', endColor, 'Marker', '^', 'LineWidth', lineWidth);
plot(xRange_Start, lambda_Start_fit, 'Color', startColor, 'LineWidth', lineWidth);
plot(xRange_Mid,   lambda_Mid_fit,   'Color', midColor,   'LineWidth', lineWidth);
plot(xRange_End,   lambda_End_fit,   'Color', endColor,   'LineWidth', lineWidth);
title({sprintf('Start: r=%.2f, p=%.2f, mean FC=%.2f',   corr_lambda_Start, corr_lambda_Start_P, nanmean(mean_Conn_Start)), ...
       sprintf('Middle: r=%.2f, p=%.2f, mean FC=%.2f', corr_lambda_Mid,   corr_lambda_Mid_P,   nanmean(mean_Conn_Mid)), ...
       sprintf('End: r=%.2f, p=%.2f, mean FC=%.2f',    corr_lambda_End,   corr_lambda_End_P,   nanmean(mean_Conn_End))});
xlabel('Mean connectivity (z)'); ylabel('Characteristic Path Length (AUC)');
set(gca, 'FontSize', 14);
% Plot correlation: Clustering Coefficient AUC vs mean connectivity
hFig4b = figure; set(hFig4b, 'color', 'w', 'Name', 'Corr_ClustCoeff');
hold on;
scatter(mean_Conn_Start, AUC_clust_coeff_HbT(designMat(:,1)), markerSize, 'MarkerEdgeColor', startColor, 'MarkerFaceColor', startColor, 'Marker', 'o', 'LineWidth', lineWidth);
scatter(mean_Conn_Mid,   AUC_clust_coeff_HbT(designMat(:,2)), markerSize, 'MarkerEdgeColor', midColor, 'Marker', 'x', 'LineWidth', lineWidth);
scatter(mean_Conn_End,   AUC_clust_coeff_HbT(designMat(:,3)), markerSize, 'MarkerEdgeColor', endColor, 'MarkerFaceColor', endColor, 'Marker', '^', 'LineWidth', lineWidth);
plot(xRange_Start, clust_Start_fit, 'Color', startColor, 'LineWidth', lineWidth);
plot(xRange_Mid,   clust_Mid_fit,   'Color', midColor,   'LineWidth', lineWidth);
plot(xRange_End,   clust_End_fit,   'Color', endColor,   'LineWidth', lineWidth);
title({sprintf('Start: r=%.2f, p=%.2e, mean FC=%.2f',   corr_clust_Start, corr_clust_Start_P, nanmean(mean_Conn_Start)), ...
       sprintf('Middle: r=%.2f, p=%.2e, mean FC=%.2f', corr_clust_Mid,   corr_clust_Mid_P,   nanmean(mean_Conn_Mid)), ...
       sprintf('End: r=%.2f, p=%.2e, mean FC=%.2f',    corr_clust_End,   corr_clust_End_P,   nanmean(mean_Conn_End))});
xlabel('Mean connectivity (z)'); ylabel('Clustering Coefficient (AUC)');
set(gca, 'FontSize', 14);
% Plot correlation: Small-World Index AUC vs mean connectivity
hFig4c = figure; set(hFig4c, 'color', 'w', 'Name', 'Corr_SWI');
hold on;
scatter(mean_Conn_Start, AUC_SWI_HbT(designMat(:,1)), markerSize, 'MarkerEdgeColor', startColor, 'MarkerFaceColor', startColor, 'Marker', 'o', 'LineWidth', lineWidth);
scatter(mean_Conn_Mid,   AUC_SWI_HbT(designMat(:,2)), markerSize, 'MarkerEdgeColor', midColor, 'Marker', 'x', 'LineWidth', lineWidth);
scatter(mean_Conn_End,   AUC_SWI_HbT(designMat(:,3)), markerSize, 'MarkerEdgeColor', endColor, 'MarkerFaceColor', endColor, 'Marker', '^', 'LineWidth', lineWidth);
plot(xRange_Start, swi_Start_fit, 'Color', startColor, 'LineWidth', lineWidth);
plot(xRange_Mid,   swi_Mid_fit,   'Color', midColor,   'LineWidth', lineWidth);
plot(xRange_End,   swi_End_fit,   'Color', endColor,   'LineWidth', lineWidth);
title({sprintf('Start: r=%.2f, p=%.2e, mean FC=%.2f',   corr_swi_Start, corr_swi_Start_P, nanmean(mean_Conn_Start)), ...
       sprintf('Middle: r=%.2f, p=%.2e, mean FC=%.2f', corr_swi_Mid,   corr_swi_Mid_P,   nanmean(mean_Conn_Mid)), ...
       sprintf('End: r=%.2f, p=%.2e, mean FC=%.2f',    corr_swi_End,   corr_swi_End_P,   nanmean(mean_Conn_End))});
xlabel('Mean connectivity (z)'); ylabel('Small-World Index (AUC)');
set(gca, 'FontSize', 14);
if printFigures
    print(hFig4a, '-dpng', fullfile('..\figures', 'corr_char_path_length.png'), '-r300');
    print(hFig4b, '-dpng', fullfile('..\figures', 'corr_clust_coeff.png'), '-r300');
    print(hFig4c, '-dpng', fullfile('..\figures', 'corr_swi.png'), '-r300');
end
% EOF
