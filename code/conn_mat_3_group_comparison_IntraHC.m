function [zMat1Mean,zMat2Mean,zMat3Mean,P_anova,P12,P13,P23] = ...
    conn_mat_3_group_comparison_IntraHC(fileGroup1,fileGroup2,fileGroup3,groupNames,saveFileName, ...
                                     chIdxL,group1Color,group2Color,group3Color,alpha)
% Computes intrahemispheric connectivity (left hemisphere) at three time-points,
% then runs repeated-measures ANOVA and Tukey–Kramer pairwise comparisons.

%% Load & gather Group 1
load(fileGroup1)  % expects zMatFDR, keepRun, runName, nRuns
z1=[];
for i=1:nRuns
    if keepRun(i)
        z1 = cat(3, z1, zMatFDR{i});
    end
end

%% Load & gather Group 2
clearvars -except z1 fileGroup2 fileGroup3 groupNames saveFileName chIdxL chIdxR group1Color group2Color group3Color alpha
load(fileGroup2)
z2=[];
for i=1:nRuns
    if keepRun(i)
        z2 = cat(3, z2, zMatFDR{i});
    end
end

%% Load & gather Group 3
clearvars -except z1 z2 fileGroup3 groupNames saveFileName chIdxL chIdxR group1Color group2Color group3Color alpha
load(fileGroup3)
z3=[];
for i=1:nRuns
    if keepRun(i)
        z3 = cat(3, z3, zMatFDR{i});
    end
end

%% Median connectivity matrices (for visualization)
zMat1Mean = nanmedian(z1,3);
zMat2Mean = nanmedian(z2,3);
zMat3Mean = nanmedian(z3,3);

%% Compute intrahemispheric connectivity per subject (mean across left-hemisphere channels)
N = size(z1,3);  % assumes same keepRun pattern
zVec1 = nan(N,1);
zVec2 = nan(N,1);
zVec3 = nan(N,1);
idx = 0;
for s = 1:N
    if keepRun(s)
        idx = idx + 1;
        % left-hemisphere submatrix
        m1 = z1(chIdxL, chIdxL, s);
        m2 = z2(chIdxL, chIdxL, s);
        m3 = z3(chIdxL, chIdxL, s);
        zVec1(idx) = nanmean(m1(:));
        zVec2(idx) = nanmean(m2(:));
        zVec3(idx) = nanmean(m3(:));
    end
end

%% Repeated-measures ANOVA
T = table(zVec1, zVec2, zVec3, 'VariableNames', groupNames);
Within = table(categorical(groupNames)','VariableNames',{'TimePoint'});
rm = fitrm(T, sprintf('%s-%s~1', groupNames{1}, groupNames{end}), 'WithinDesign', Within);
R = ranova(rm);
P_anova = R.pValue(1);

%% Pairwise Tukey–Kramer comparisons
C = multcompare(rm, 'TimePoint', 'ComparisonType','tukey-kramer','Alpha',alpha);
P12 = nan; P13 = nan; P23 = nan;
for k = 1:height(C)
    t1 = string(C.TimePoint_1(k)); t2 = string(C.TimePoint_2(k)); p = C.pValue(k);
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
    'zMat1Mean','zMat2Mean','zMat3Mean', ...
    'P_anova','P12','P13','P23');

%% Boxplot of intra-HC across sessions
allData = {zVec1; zVec2; zVec3};
groupVec = [ones(idx,1); 2*ones(idx,1); 3*ones(idx,1)];
colors = {group1Color, group2Color, group3Color};
markers = {'o','x','s'};

hFig = figure; set(hFig,'Color','w');
boxplot(cell2mat(allData), groupVec,'Colors',[0 0 0]); hold on;
for g = 1:3
    scatter(repmat(g,numel(allData{g}),1), allData{g}, markers{g}, ...
            'MarkerEdgeColor', colors{g}, 'jitter','on','jitterAmount',0.15);
end
set(gca,'XTick',1:3,'XTickLabel',groupNames,'FontSize',14);
ylabel('IntraHC (z-score)');
title(sprintf('RM-ANOVA p=%.3f | 1v2 p=%.3f, 1v3 p=%.3f, 2v3 p=%.3f', ...
    P_anova, P12, P13, P23));

%% Save figure
set(hFig,'Units','inches','Position',[0.1 0.1 6 6],'PaperPosition',[0. 0.1 6 6]);
print(hFig,'-dpng',fullfile('..','figures',[saveFileName '_IntraHC.png']),'-r300');
fprintf('Intrahemispheric RM-ANOVA done (p=%.3f)\n',P_anova);
end
