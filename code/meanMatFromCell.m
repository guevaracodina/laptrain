function meanMat = meanMatFromCell(zMatCell, channelIdx, alpha)
% meanMatFromCell   Extracts sub‐matrices, tests each entry vs. zero, and
%                  returns only those means significantly ≠0.
% 
% Inputs:
%   zMatCell  – cell array of full matrices
%   channelIdx– vector of channel indices to keep
%   alpha     – significance threshold (default 0.05)
% Output:
%   meanMat   – matrix of mean values, with non‐significant entries set to NaN

    if nargin<3, alpha = 0.05; end

    % 1) extract the sub‐matrices
    subMats = cellfun(@(X) X(channelIdx,channelIdx), zMatCell, 'UniformOutput',false);

    % 2) stack into a 3‐D array: rows×cols×runs
    stack  = cat(3, subMats{:});

    % 3) one‐sample t‐test vs. zero along the 3rd dim
    [h, ~] = ttest(stack, 0, 'Alpha', alpha, 'Dim', 3);

    % 4) logical mask of “significantly non‐zero” (h==1 yields false for NaNs)
    sigMask = (h == 1);

    % 5) mean across runs (ignoring NaNs)
    meanAll = nanmean(stack,3);

    % 6) mask out non‐significant (and NaN) entries
    meanMat = meanAll;
    meanMat(~sigMask) = NaN;
end


% ONLY AVERAGE
% function meanMat = meanMatFromCell(zMatCell, channelIdx)
% % 1) extract the 22×22 sub‐matrices
% subMats = cellfun(@(X) X(channelIdx,channelIdx), zMatCell, 'UniformOutput',false);
% % 2) stack into a 3D array
% stack  = cat(3, subMats{:});
% % 3) mean‐across runs (ignoring NaNs)
% meanMat = nanmean(stack,3);
% % 4) overall scalar mean
% % meanValue = nanmean(meanMat(:));
% end