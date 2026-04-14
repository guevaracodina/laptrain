function [zMat1Mean,zMat2Mean,zMat3Mean,rMat1Mean,rMat2Mean,rMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_InterHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                       chIdxL,chIdxR,group1Color,group2Color,group3Color,alpha)
% Computes inter-hemispheric connectivity at three time-points (same subjects),
% then runs repeated-measures ANOVA and Tukey–Kramer pairwise comparisons.

%% Load & gather Group 1
load(fileGroup1)  % expects zMatFDR, rMatFDR, keepRun, runName, nRuns
z1=[]; r1=[]; names1={};
for i=1:nRuns
    if keepRun(i)
        z1=cat(3,z1,zMatFDR{i});
        r1=cat(3,r1,zMatFDR{i});
        names1{end+1,1}=runName{i};
    end
end

%% Load & gather Group 2
clearvars -except z1 r1 names1 fileGroup2 fileGroup3 groupNames saveFileName chIdxL chIdxR group1Color group2Color group3Color alpha
load(fileGroup2)
z2=[]; r2=[]; names2={};
for i=1:nRuns
    if keepRun(i)
        z2=cat(3,z2,zMatFDR{i});
        r2=cat(3,r2,zMatFDR{i});
        names2{end+1,1}=runName{i};
    end
end

%% Load & gather Group 3
clearvars -except z1 r1 names1 z2 r2 names2 fileGroup3 groupNames saveFileName chIdxL chIdxR group1Color group2Color group3Color alpha
load(fileGroup3)
z3=[]; r3=[]; names3={};
for i=1:nRuns
    if keepRun(i)
        z3=cat(3,z3,zMatFDR{i});
        r3=cat(3,r3,zMatFDR{i});
        names3{end+1,1}=runName{i};
    end
end

%% Median connectivity matrices
zMat1Mean = nanmedian(z1,3);
zMat2Mean = nanmedian(z2,3);
zMat3Mean = nanmedian(z3,3);
rMat1Mean = nanmedian(r1,3);
rMat2Mean = nanmedian(r2,3);
rMat3Mean = nanmedian(r3,3);

%% Extract inter-HC per subject
nCh = numel(chIdxL);
% Preallocate subject connectivity vectors
zVec1 = nan(size(z1,3), nCh);
zVec2 = nan(size(z2,3), nCh);
zVec3 = nan(size(z3,3), nCh);

for s=1:size(z1,3)
    for c=1:nCh
        zVec1(s,c) = z1(chIdxL(c), chIdxR(c), s);
    end
end
for s=1:size(z2,3)
    for c=1:nCh
        zVec2(s,c) = z2(chIdxL(c), chIdxR(c), s);
    end
end
for s=1:size(z3,3)
    for c=1:nCh
        zVec3(s,c) = z3(chIdxL(c), chIdxR(c), s);
    end
end

%% Median per subject
zSubj1 = nanmedian(zVec1,2);
zSubj2 = nanmedian(zVec2,2);
zSubj3 = nanmedian(zVec3,2);

%% Match subject names across time-points
% [common12,i1,i2] = intersect(names1, names2);
% [commonAll,i12,i3] = intersect(common12, names3);
% i1 = i1(i12);  i2 = i2(i12);  i3 = i3;
% zSubj1 = zSubj1(i1);
% zSubj2 = zSubj2(i2);
% zSubj3 = zSubj3(i3);

%% Build table for RM ANOVA
T = table(zSubj1, zSubj2, zSubj3, 'VariableNames', groupNames);
Within = table(categorical(groupNames)', 'VariableNames', {'TimePoint'});

%% Fit RM model
rm = fitrm(T, sprintf('%s-%s~1', groupNames{1}, groupNames{end}), 'WithinDesign', Within);
R = ranova(rm);
P_anova = R.pValue(1);

%% Pairwise Tukey–Kramer
C = multcompare(rm, 'TimePoint', 'ComparisonType', 'tukey-kramer', 'Alpha', alpha);
P12 = nan; P13 = nan; P23 = nan;
for k = 1:height(C)
    t1 = string(C.TimePoint_1(k));
    t2 = string(C.TimePoint_2(k));
    p  = C.pValue(k);
    if t1 == groupNames{1} && t2 == groupNames{2}
        P12 = p;
    elseif t1 == groupNames{1} && t2 == groupNames{3}
        P13 = p;
    elseif t1 == groupNames{2} && t2 == groupNames{3}
        P23 = p;
    end
end

%% Save results
save(fullfile('..','data',[saveFileName '.mat']), ...
    'zMat1Mean', 'zMat2Mean', 'zMat3Mean', ...
    'rMat1Mean', 'rMat2Mean', 'rMat3Mean', ...
    'P_anova', 'P12', 'P13', 'P23');

%% Boxplot
allData = {zSubj1; zSubj2; zSubj3};
groupVec = [ones(size(zSubj1)); 2*ones(size(zSubj2)); 3*ones(size(zSubj3))];
colors = {group1Color, group2Color, group3Color};
markers = {'o', 'x', 's'};

hFig = figure; set(hFig, 'Color', 'w', 'Name', 'Inter-HC');
boxplot(cell2mat(allData), groupVec, 'Colors', [0 0 0]); hold on;
for g = 1:3
    scatter(repmat(g, numel(allData{g}), 1), allData{g}, markers{g}, ...
            'MarkerEdgeColor', colors{g}, 'jitter', 'on', 'jitterAmount', 0.15);
end
set(gca, 'XTick', 1:3, 'XTickLabel', groupNames, 'FontSize', 14);
ylabel('InterHC (z-score)');
title(sprintf('RM-ANOVA p=%.3f | 1v2 p=%.3f, 1v3 p=%.3f, 2v3 p=%.3f', ...
    P_anova, P12, P13, P23));

%% Save figure
set(hFig, 'Units', 'inches', 'Position', [0.1 0.1 6 6], 'PaperPosition', [0.1 0.1 6 6]);
print(hFig, '-dpng', fullfile('..','figures',[saveFileName '_InterHC.png']), '-r300');
fprintf('Inter-hemispherical RM-ANOVA done (p=%.3f)\n', P_anova);
end
