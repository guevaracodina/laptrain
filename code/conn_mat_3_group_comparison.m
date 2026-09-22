function [zMat1Mean,zMat2Mean,zMat3Mean, ...
          P12,P13,P23,pFDR12,pFDR13,pFDR23, ...
          T12,T13,T23] = ...
    conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
                                groupNames, saveFileName, alpha, chIdxL, chIdxR, symmColorZero)
% 3-group RM-ANOVA + pairwise p- and t-value matrices, 4-row plotting.
% If symmColorZero is true, rows 1 (mean Z) and 4 (t-statistic) use
% symmetric color limits around zero.
fprintf('Computing 3-group (%s %s %s) RM-ANOVA\n + pairwise p- and t-value matrices\n', ...
        fileGroup1, fileGroup2, fileGroup3)

if nargin<9
    symmColorZero = false;
end

%% Load & gather each group
for g = 1:3
    file = eval(sprintf('fileGroup%d',g));
    clear zMatFDR keepRun runName nRuns
    load(file)  % loads zMatFDR, keepRun, runName, nRuns
    z = [];
    for i=1:nRuns
        if keepRun(i)
            z = cat(3, z, zMatFDR{i});
        end
    end
    eval(sprintf('zMat%d = z;',g));
end

%% Means
zMat1Mean = nanmean(zMat1,3);
zMat2Mean = nanmean(zMat2,3);
zMat3Mean = nanmean(zMat3,3);

%% Prepare stats matrices & basic limits
nCh = size(zMat1,1);
P12 = nan(nCh); P13 = nan(nCh); P23 = nan(nCh);
T12 = nan(nCh); T13 = nan(nCh); T23 = nan(nCh);

minVal = min([zMat1Mean(:); zMat2Mean(:); zMat3Mean(:)]);
maxVal = max([zMat1Mean(:); zMat2Mean(:); zMat3Mean(:)]);


%% Color limits for means (row1)
if symmColorZero
    limMean = max(abs([minVal, maxVal]));
    clim_mean = [-limMean, limMean];
else
    clim_mean = [minVal, maxVal];
end

%% Stats per edge
for iCh=1:nCh
  for jCh=1:nCh
    x1 = squeeze(zMat1(iCh,jCh,:)); x2 = squeeze(zMat2(iCh,jCh,:)); x3 = squeeze(zMat3(iCh,jCh,:));
    if any([all(isnan(x1)), all(isnan(x2)), all(isnan(x3))]) || ...
       ~(numel(x1)==numel(x2) && numel(x2)==numel(x3))
       continue;
    end
    % ANOVA & p-values
    T = table(x1,x2,x3,'VariableNames',groupNames);
    W = table(categorical(groupNames)','VariableNames',{'Cond'});
    rm = fitrm(T, sprintf('%s-%s~1',groupNames{1},groupNames{end}),'WithinDesign',W);
    C = multcompare(rm,'Cond','ComparisonType','tukey-kramer','Alpha',alpha);
    for k=1:height(C)
      n1=string(C.Cond_1(k)); n2=string(C.Cond_2(k)); p=C.pValue(k);
      if n1==groupNames{1}&&n2==groupNames{2}, P12(iCh,jCh)=p;
      elseif n1==groupNames{1}&&n2==groupNames{3}, P13(iCh,jCh)=p;
      elseif n1==groupNames{2}&&n2==groupNames{3}, P23(iCh,jCh)=p; end
    end
    % pairwise t-tests
    [~,~,~,s12] = ttest(x1,x2); T12(iCh,jCh)=s12.tstat;
    [~,~,~,s13] = ttest(x1,x3); T13(iCh,jCh)=s13.tstat;
    [~,~,~,s23] = ttest(x2,x3); T23(iCh,jCh)=s23.tstat;
  end
end

%% FDR corrections
[pFDR12,~] = conn_mat_fdr(P12,alpha);
[pFDR13,~] = conn_mat_fdr(P13,alpha);
[pFDR23,~] = conn_mat_fdr(P23,alpha);

%% Prepare color limits for t-statistics (row4)
tmin = min([T12(:); T13(:); T23(:)]);
tmax = max([T12(:); T13(:); T23(:)]);
if symmColorZero
    limT = max(abs([tmin, tmax]));
    clim_t = [-limT, limT];
else
    clim_t = [tmin, tmax];
end

%% Save results
save(fullfile('..','data',[saveFileName '.mat']), ...
    'zMat1Mean','zMat2Mean','zMat3Mean', ...
    'P12','P13','P23','pFDR12','pFDR13','pFDR23', ...
    'T12','T13','T23');

%% Plot 4x3
chIdx = reshape([chIdxL(:)'; chIdxR(:)'],[],1);
axisLabels = [chIdxL chIdxR];

MNI_brain_regions = {...
    {'Frontal_Inf_Tri_R'   }
    {'Frontal_Inf_Tri_R'   }
    {'SSC'}
    {'Frontal_Inf_Tri_R'   }
    {'Frontal_Sup_R'       }
    {'Frontal_Mid_R'       }
    {'Frontal_Mid_R'       }
    {'Frontal_Sup_R'       }
    {'Frontal_Sup_R'       }
    {'Frontal_Sup_R'       }
    {'Frontal_Sup_R'       }
    {'Frontal_Sup_Medial_L'}
    {'Frontal_Sup_Medial_L'}
    {'SSC'}
    {'Frontal_Sup_R'       }
    {'Frontal_Sup_L'       }
    {'Frontal_Mid_L'       }
    {'Frontal_Sup_L'       }
    {'Frontal_Mid_L'       }
    {'Frontal_Inf_Tri_L'   }
    {'Frontal_Sup_Mid_L'   }
    {'Frontal_Mid_L'       }
    {'Frontal_Inf_Tri_L'   }
    {'Frontal_Inf_Tri_L'   }};

regionLabels = MNI_brain_regions([chIdxL'; chIdxR']);
labels = cellfun(@(c)c{1}, regionLabels, 'UniformOutput', false);
% p-values matrix color limits
thr12 = -log10(P12); thr12(P12>=alpha)=nan;
thr13 = -log10(P13); thr13(P13>=alpha)=nan;
thr23 = -log10(P23); thr23(P23>=alpha)=nan;
clim_p = [min([thr12(:);thr13(:);thr23(:)],[],'omitnan'), max([thr12(:);thr13(:);thr23(:)],[],'omitnan')];

hFig = figure('Color','w');

% Row 1: mean Z-matrices
mats = {zMat1Mean, zMat2Mean, zMat3Mean};
Tmats = {T12,T13,T23};

for k=1:3
  ax = subplot(4,3,k);
  imagesc(mats{k}(chIdx,chIdx), clim_mean);
  title(groupNames{k});
  axis image;
  c = colorbar; c.Label.String = 'z-score';
  colormap(ax, ioi_get_colormap('redbluecmap'));
  set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end), ...
         'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
end

% row2 uncorrected p
titles_unc = { ...
  sprintf('%s vs %s', groupNames{1}, groupNames{2}), ...
  sprintf('%s vs %s', groupNames{1}, groupNames{3}), ...
  sprintf('%s vs %s', groupNames{2}, groupNames{3}) ...
};
for k=1:3
  ax=subplot(4,3,3+k);
  mat = {thr12,thr13,thr23}; m = mat{k};
  imagesc(m(chIdx,chIdx),clim_p); title(titles_unc{k});
  axis image; a = colorbar; a.Label.String = '-log_{10}(p)'; colormap(ax,flipud(ioi_get_colormap('linearl')));
  set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
         'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
end
% row3 FDR
titles_fdr = strrep(titles_unc,'',' FDR');
for k=1:3
  ax=subplot(4,3,6+k);
  mat = {pFDR12,pFDR13,pFDR23}; m = mat{k}; m(m>=alpha)=nan;
  imagesc(-log10(m(chIdx,chIdx)),clim_p); title(titles_fdr{k});
  axis image; a = colorbar; a.Label.String = '-log_{10}(q)'; colormap(ax,flipud(ioi_get_colormap('linearl')));
  set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
         'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
end

% Row 4: t-statistics
for k=1:3
  ax = subplot(4,3,9+k);
  imagesc(Tmats{k}(chIdx,chIdx), clim_t);
  title(strrep(titles_unc{k},'',' t-stat'));
  axis image;
  c = colorbar; c.Label.String = 't-statistic';
  colormap(ax, ioi_get_colormap('redbluecmap'));
  set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end), ...
         'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
  % Brain regions are not clearly visible, hence their labels were removed
  % if k==1
  %   set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'YTick',1:numel(labels),'YTickLabel',labels,'FontSize',6);
  %   xtickangle(ax,45);
  % end
end

%% Finalize and save figure
set(hFig,'Units','inches','Position',[0.1 0.1 9 12],'PaperPosition',[0.1 0.1 9 12]);
print(hFig,'-dpng',fullfile('..','figures',[saveFileName '_conn_mat4row.png']),'-r300');
disp('Finished 3-group comparison!')
end

% WORKING FUNCTION! no symmetrical colormaps
% function [zMat1Mean,zMat2Mean,zMat3Mean, ...
%           P12,P13,P23,pFDR12,pFDR13,pFDR23, ...
%           T12,T13,T23] = ...
%     conn_mat_3_group_comparison(fileGroup1, fileGroup2, fileGroup3, ...
%                                 groupNames, saveFileName, alpha, chIdxL, chIdxR)
% % 3-group RM-ANOVA + pairwise p- and t-value matrices, 4-row plotting.
% fprintf('Computing 3-group (%s %s %s) RM-ANOVA\n + pairwise p- and t-value matrices\n', fileGroup1, fileGroup2, fileGroup3)
% %% Load & gather each group
% for g = 1:3
%     file = eval(sprintf('fileGroup%d',g));
%     clear zMat rMat names keepRun runName nRuns iHb
%     load(file)  % loads zMatFDR, keepRun, runName, nRuns
%     z = [];
%     for i=1:nRuns
%         if keepRun(i)
%             z = cat(3,z, zMatFDR{i});
%         end
%     end
%     eval(sprintf('zMat%d = z;',g));
% end
% 
% %% Means
% zMat1Mean = nanmean(zMat1,3);
% zMat2Mean = nanmean(zMat2,3);
% zMat3Mean = nanmean(zMat3,3);
% 
% %% Prepare
% nCh = size(zMat1,1);
% P12 = nan(nCh); P13 = nan(nCh); P23 = nan(nCh);
% T12 = nan(nCh); T13 = nan(nCh); T23 = nan(nCh);
% minVal = min([zMat1Mean(:);zMat2Mean(:);zMat3Mean(:)]);
% maxVal = max([zMat1Mean(:);zMat2Mean(:);zMat3Mean(:)]);
% 
% %% Stats per edge
% for iCh=1:nCh
%   for jCh=1:nCh
%     x1 = squeeze(zMat1(iCh,jCh,:)); x2 = squeeze(zMat2(iCh,jCh,:)); x3 = squeeze(zMat3(iCh,jCh,:));
%     if any([ all(isnan(x1)), all(isnan(x2)), all(isnan(x3)) ]) || ...
%        ~(numel(x1)==numel(x2) && numel(x2)==numel(x3))
%        continue; end
%     % ANOVA & p-values
%     T = table(x1,x2,x3,'VariableNames',groupNames);
%     W = table(categorical(groupNames)','VariableNames',{'Cond'});
%     rm = fitrm(T, sprintf('%s-%s~1',groupNames{1},groupNames{end}),'WithinDesign',W);
%     C = multcompare(rm,'Cond','ComparisonType','tukey-kramer','Alpha',alpha);
%     for k=1:height(C)
%       n1=string(C.Cond_1(k)); n2=string(C.Cond_2(k)); p=C.pValue(k);
%       if n1==groupNames{1}&&n2==groupNames{2}, P12(iCh,jCh)=p;
%       elseif n1==groupNames{1}&&n2==groupNames{3}, P13(iCh,jCh)=p;
%       elseif n1==groupNames{2}&&n2==groupNames{3}, P23(iCh,jCh)=p; end
%     end
%     % pairwise t-tests
%     [~,~,~,s12] = ttest(x1,x2); T12(iCh,jCh)=s12.tstat;
%     [~,~,~,s13] = ttest(x1,x3); T13(iCh,jCh)=s13.tstat;
%     [~,~,~,s23] = ttest(x2,x3); T23(iCh,jCh)=s23.tstat;
%   end
% end
% 
% %% FDR
% [pFDR12,~]=conn_mat_fdr(P12,alpha);
% [pFDR13,~]=conn_mat_fdr(P13,alpha);
% [pFDR23,~]=conn_mat_fdr(P23,alpha);
% 
% %% Save
% save(fullfile('..','data',[saveFileName '.mat']), ...
%     'zMat1Mean','zMat2Mean','zMat3Mean', ...
%     'P12','P13','P23','pFDR12','pFDR13','pFDR23', ...
%     'T12','T13','T23');
% 
% %% Plot 4x3
% chIdx = reshape([chIdxL(:)';chIdxR(:)'],[],1);
% axisLabels = [chIdxL chIdxR];
% MNI_brain_regions = {...
%     {'Frontal_Inf_Tri_R'   }
%     {'Frontal_Inf_Tri_R'   }
%     {'SSC'}
%     {'Frontal_Inf_Tri_R'   }
%     {'Frontal_Sup_R'       }
%     {'Frontal_Mid_R'       }
%     {'Frontal_Mid_R'       }
%     {'Frontal_Sup_R'       }
%     {'Frontal_Sup_R'       }
%     {'Frontal_Sup_R'       }
%     {'Frontal_Sup_R'       }
%     {'Frontal_Sup_Medial_L'}
%     {'Frontal_Sup_Medial_L'}
%     {'SSC'}
%     {'Frontal_Sup_R'       }
%     {'Frontal_Sup_L'       }
%     {'Frontal_Mid_L'       }
%     {'Frontal_Sup_L'       }
%     {'Frontal_Mid_L'       }
%     {'Frontal_Inf_Tri_L'   }
%     {'Frontal_Sup_Mid_L'   }
%     {'Frontal_Mid_L'       }
%     {'Frontal_Inf_Tri_L'   }
%     {'Frontal_Inf_Tri_L'   }};
% 
% regionLabels = MNI_brain_regions([chIdxL'; chIdxR']);
% 
% % color limits
% thr12 = -log10(P12); thr12(P12>=alpha)=nan;
% thr13 = -log10(P13); thr13(P13>=alpha)=nan;
% thr23 = -log10(P23); thr23(P23>=alpha)=nan;
% clim_p = [min([thr12(:);thr13(:);thr23(:)],[],'omitnan'), max([thr12(:);thr13(:);thr23(:)],[],'omitnan')];
% clim_t = [min([T12(:);T13(:);T23(:)]), max([T12(:);T13(:);T23(:)])];
% 
% hFig=figure('Color','w');
% % row1 means
% mats={zMat1Mean,zMat2Mean,zMat3Mean};
% Tmats = {T12,T13,T23};
% for k=1:3
%   ax=subplot(4,3,k);
%   imagesc(mats{k}(chIdx,chIdx),[minVal maxVal]); title(groupNames{k});
%   axis image; a = colorbar; a.Label.String = 'z-score'; colormap(ax,ioi_get_colormap('redbluecmap'));
%   set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
%          'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
% end
% % row2 uncorrected p
% titles_unc = { ...
%   sprintf('%s vs %s', groupNames{1}, groupNames{2}), ...
%   sprintf('%s vs %s', groupNames{1}, groupNames{3}), ...
%   sprintf('%s vs %s', groupNames{2}, groupNames{3}) ...
% };
% for k=1:3
%   ax=subplot(4,3,3+k);
%   mat = {thr12,thr13,thr23}; m = mat{k};
%   imagesc(m(chIdx,chIdx),clim_p); title(titles_unc{k});
%   axis image; a = colorbar; a.Label.String = '-log_{10}(p)'; colormap(ax,flipud(ioi_get_colormap('linearl')));
%   set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
%          'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
% end
% % row3 FDR
% titles_fdr = strrep(titles_unc,'',' FDR');
% for k=1:3
%   ax=subplot(4,3,6+k);
%   mat = {pFDR12,pFDR13,pFDR23}; m = mat{k}; m(m>=alpha)=nan;
%   imagesc(-log10(m(chIdx,chIdx)),clim_p); title(titles_fdr{k});
%   axis image; a = colorbar; a.Label.String = '-log_{10}(q)'; colormap(ax,flipud(ioi_get_colormap('linearl')));
%   set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
%          'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
% end
% % row4 t-values
% titles_t = strrep(titles_unc,'',' t');
% labels = cellfun(@(c)c{1}, MNI_brain_regions, 'UniformOutput', false);
% 
% for k=1:3
%   ax = subplot(4,3,9+k);
%   imagesc(Tmats{k}(chIdx,chIdx),clim_t);
%   title(titles_t{k});
% 
%   axis image; 
%   a = colorbar; 
%   a.Label.String = 't-statistic'; 
%   colormap(ax,ioi_get_colormap('redbluecmap'));
% 
%   set(ax,'XTick',1:2:numel(axisLabels),'XTickLabel',axisLabels(1:2:end),...
%          'YTick',1:2:numel(axisLabels),'YTickLabel',axisLabels(1:2:end));
% 
%   if k == 1
%   % use the MNI region names every other tick
%   set(ax, ...
%       'XTick',    1:numel(labels),         'XTickLabel', labels(1:end), ...
%       'YTick',    1:numel(labels),         'YTickLabel', labels(1:end));
%   set(ax, 'FontSize', 6);
%   xtickangle(ax,45);
%   end
% end
% 
% %% Save
% set(hFig,'Units','inches','Position',[0.1 0.1 9 12],'PaperPosition',[0.1 0.1 9 12]);
% print(hFig,'-dpng',fullfile('..','figures',[saveFileName '_conn_mat4row.png']),'-r300');
% disp('Finished 3-group comparison!')
% end

    
    