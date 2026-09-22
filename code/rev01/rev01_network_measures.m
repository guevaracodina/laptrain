function M = rev01_network_measures(dataDir, condition, chromophore, nPerms, seed)
% REV01_NETWORK_MEASURES  Graph metrics and repeated-measures tests, revision 01.
%
%   M = rev01_network_measures(dataDir, condition, chromophore, nPerms, seed)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Reproduces the network analysis of fc_network_measures_script_2025.m with
%   two changes demanded by the review:
%     (1) it runs for any chromophore, so that HbO, HbR and HbT can be
%         compared (Associate Editor major concern 1b);
%     (2) the group tests respect the repeated-measures design, using
%         permutationTestRM3 and permutationTestRMpair rather than the
%         independent-sample permutationTest3/permutationTest
%         (Associate Editor major concern 1a);
%     (3) the clustering coefficient is computed on a BINARIZED graph, as
%         clustering_coef_bu requires and as the manuscript states. The
%         submitted pipeline passed the weighted output of
%         threshold_proportional, so sum(S(:))/(k^2-k) summed correlation
%         weights rather than counting triangles and the metric scaled
%         linearly with connectivity strength. Both versions are returned:
%         'clust' is the corrected binary clustering coefficient and
%         'clustW' reproduces the submitted quantity, so the difference
%         between them can be inspected directly.
%     Characteristic path length is likewise returned both ways, weighted
%     ('lambda', as submitted) and binary ('lambdaB'). The small-world index
%     is left exactly as submitted, computed by swi.m on the weighted graph.
%   It also returns the per-subject metrics so that data-quality covariates
%   can be added downstream (concern 1d).
%
%   INPUT
%       dataDir     path to the project 'data' folder, e.g. '..\data'
%       condition   'Lap' or 'Resting'
%       chromophore 'HbO', 'HbR' or 'HbT'
%       nPerms      permutations for the group tests (default 5000)
%       seed        rng seed for reproducibility (default 42)
%
%   OUTPUT struct M
%       threshold   [1 x nThr] sparsity thresholds
%       lambda      [nSub x nThr x 3] characteristic path length
%       clust       [nSub x nThr x 3] mean clustering coefficient
%       SWI         [nSub x nThr x 3] small-world index
%       AUClambda, AUCclust, AUCswi   [nSub x 3] area under the curve
%       meanFC      [nSub x 3] mean of the off-diagonal Fisher z values
%       valid       [nSub x 3] logical, run retained by conn_mat_subject_level
%       pAUC*       global repeated-measures permutation p-values
%       pAUC*_12/13/23   post-hoc paired permutation p-values
%       pThr*       [1 x nThr] threshold-wise global p-values
%
%   Requires the Brain Connectivity Toolbox on the path.

if nargin < 4 || isempty(nPerms), nPerms = 5000; end
if nargin < 5 || isempty(seed),   seed   = 42;   end
rng(seed);

tags = {'begin', 'mid', 'end'};
M.condition   = condition;
M.chromophore = chromophore;

% Long channels of the prefrontal montage; short channels 3 and 24 excluded
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];
nNode = numel(channelIdx);

%% Load the three time points, keeping subjects aligned across time
Z = [];  valid = [];
for iTime = 1:3
    f = fullfile(dataDir, sprintf('%s%sConn%s.mat', tags{iTime}, condition, chromophore));
    if ~isfile(f)
        error('rev01_network_measures:missingFile', 'Not found: %s', f);
    end
    S = load(f, 'zMatFDR', 'keepRun');
    nSub = numel(S.zMatFDR);
    if isempty(Z)
        Z = nan(nNode, nNode, nSub, 3);
        valid = false(nSub, 3);
    end
    for iSub = 1:nSub
        z = S.zMatFDR{iSub};
        if isempty(z) || ~isequal(size(z), [24 24])
            continue                     % run dropped upstream; stays NaN
        end
        Z(:, :, iSub, iTime) = z(channelIdx, channelIdx);
        valid(iSub, iTime)   = S.keepRun(iSub);
    end
end
[~, ~, nSub, ~] = size(Z);
M.valid = valid;

%% Mean connectivity per subject and time point (off-diagonal only)
offDiag = ~eye(nNode);
meanFC  = nan(nSub, 3);
for iTime = 1:3
    for iSub = 1:nSub
        z = Z(:, :, iSub, iTime);
        meanFC(iSub, iTime) = mean(z(offDiag), 'omitnan');
    end
end
M.meanFC = meanFC;

%% Graph metrics across the sparsity range used in the submitted analysis
threshold = 0.1:0.01:0.34;
nThr      = numel(threshold);
M.threshold = threshold;

lambda  = nan(nSub, nThr, 3);   % weighted characteristic path length, as submitted
lambdaB = nan(nSub, nThr, 3);   % binary characteristic path length
clust   = nan(nSub, nThr, 3);   % CORRECTED: binary clustering coefficient
clustW  = nan(nSub, nThr, 3);   % as submitted: weight sum, kept for comparison
SWIm    = nan(nSub, nThr, 3);   % small-world index, swi.m unchanged

for iTime = 1:3
    for iThr = 1:nThr
        for iSub = 1:nSub
            W = Z(:, :, iSub, iTime);
            if all(isnan(W(:))), continue; end
            Wt = threshold_proportional(W, threshold(iThr));
            Wb = weight_conversion(Wt, 'binarize');

            % Weighted path length, as in the submitted analysis
            D = distance_wei(weight_conversion(Wt, 'lengths'));
            lambda(iSub, iThr, iTime) = charpath(D, 0, 0);

            % Binary path length, consistent with a binarized network
            DB = distance_wei(weight_conversion(Wb, 'lengths'));
            lambdaB(iSub, iThr, iTime) = charpath(DB, 0, 0);

            % clustering_coef_bu requires a binary matrix
            clust(iSub, iThr, iTime)  = mean(clustering_coef_bu(Wb), 'omitnan');
            clustW(iSub, iThr, iTime) = mean(clustering_coef_bu(Wt), 'omitnan');

            % Unchanged from the submitted analysis
            SWIm(iSub, iThr, iTime) = swi(Wt);
        end
    end
    fprintf('  %s %s: time point %d of 3 done\n', condition, chromophore, iTime);
end
M.lambda = lambda;  M.lambdaB = lambdaB;
M.clust  = clust;   M.clustW  = clustW;
M.SWI    = SWIm;

%% Area under the curve over the sparsity range
% At low sparsity the thresholded graph can fragment, in which case charpath
% and therefore the small-world index are undefined. A plain trapz would turn
% the whole AUC into NaN for that subject, which silently drops subjects from
% the group test. The AUC is integrated over the thresholds at which the
% metric is defined, and the coverage is recorded so that it can be reported.
[M.AUClambda,  M.covLambda]  = aucOmitNaN(threshold, lambda);
[M.AUClambdaB, M.covLambdaB] = aucOmitNaN(threshold, lambdaB);
[M.AUCclust,   M.covClust]   = aucOmitNaN(threshold, clust);
[M.AUCclustW,  M.covClustW]  = aucOmitNaN(threshold, clustW);
[M.AUCswi,     M.covSWI]     = aucOmitNaN(threshold, SWIm);

%% Repeated-measures permutation tests on the AUC
[M.pAUClambda,  M.dAUClambda,  M.FAUClambda]  = permutationTestRM3(M.AUClambda,  nPerms);
[M.pAUClambdaB, M.dAUClambdaB, M.FAUClambdaB] = permutationTestRM3(M.AUClambdaB, nPerms);
[M.pAUCclust,   M.dAUCclust,   M.FAUCclust]   = permutationTestRM3(M.AUCclust,   nPerms);
[M.pAUCclustW,  M.dAUCclustW,  M.FAUCclustW]  = permutationTestRM3(M.AUCclustW,  nPerms);
[M.pAUCswi,     M.dAUCswi,     M.FAUCswi]     = permutationTestRM3(M.AUCswi,     nPerms);

pairs = {[1 2], [1 3], [2 3]};
names = {'12', '13', '23'};
for k = 1:3
    a = pairs{k}(1); b = pairs{k}(2);
    for mv = {'lambda', 'lambdaB', 'clust', 'clustW', 'swi'}
        fld = ['AUC' mv{1}];
        if strcmp(mv{1}, 'swi'), fld = 'AUCswi'; end
        [M.(['pAUC' mv{1} '_' names{k}]), M.(['dAUC' mv{1} '_' names{k}])] = ...
            permutationTestRMpair(M.(fld)(:, a), M.(fld)(:, b), nPerms);
    end
end

%% Threshold-wise global tests, for the starred markers in Figure 2a
M.pThrLambda = nan(1, nThr);
M.pThrClust  = nan(1, nThr);
M.pThrClustW = nan(1, nThr);
M.pThrSWI    = nan(1, nThr);
for iThr = 1:nThr
    M.pThrLambda(iThr) = permutationTestRM3(squeeze(lambda(:, iThr, :)), nPerms);
    M.pThrClust(iThr)  = permutationTestRM3(squeeze(clust(:,  iThr, :)), nPerms);
    M.pThrClustW(iThr) = permutationTestRM3(squeeze(clustW(:, iThr, :)), nPerms);
    M.pThrSWI(iThr)    = permutationTestRM3(squeeze(SWIm(:,   iThr, :)), nPerms);
end

%% Association between mean connectivity and network topology, per time point
for iTime = 1:3
    ok = ~isnan(meanFC(:, iTime)) & ~isnan(M.AUCclust(:, iTime));
    M.rFCclust(iTime)  = corr(meanFC(ok, iTime), M.AUCclust(ok, iTime));
    okW = ~isnan(meanFC(:, iTime)) & ~isnan(M.AUCclustW(:, iTime));
    M.rFCclustW(iTime) = corr(meanFC(okW, iTime), M.AUCclustW(okW, iTime));
    M.rFClambda(iTime) = corr(meanFC(ok, iTime), M.AUClambda(ok, iTime));
    M.rFCswi(iTime)    = corr(meanFC(ok, iTime), M.AUCswi(ok, iTime));
end

end

%% ------------------------------------------------------------------------
function [A, cov] = aucOmitNaN(thr, X)
% Trapezoidal AUC over the thresholds where X is finite.
% X is [nSub x nThr x nTime]; A and cov are [nSub x nTime].
[nSub, ~, nTime] = size(X);
A   = nan(nSub, nTime);
cov = zeros(nSub, nTime);
span = thr(end) - thr(1);
for iTime = 1:nTime
    for iSub = 1:nSub
        v  = squeeze(X(iSub, :, iTime));
        ok = isfinite(v);
        cov(iSub, iTime) = sum(ok) / numel(v);
        if sum(ok) < 2, continue; end
        % Rescale to the full threshold span so that subjects with different
        % coverage remain on a comparable scale.
        a = trapz(thr(ok), v(ok));
        A(iSub, iTime) = a * span / (thr(find(ok, 1, 'last')) - thr(find(ok, 1)));
    end
end
end

% EOF
