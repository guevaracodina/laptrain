function I = rev01_threshold_invariance(dataDir, condition, chromophore, nSubUse, seed)
% REV01_THRESHOLD_INVARIANCE  Is clustering sensitive to global FC inflation?
%
%   I = rev01_threshold_invariance(dataDir, condition, chromophore, nSubUse, seed)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Associate Editor major concern 1d argues that poor scalp coupling and
%   residual motion inflate correlation, and that this could inflate
%   density-dependent clustering, so that the mid-shift clustering effect
%   might be a data-quality artefact.
%
%   The network metrics here are computed on graphs built with
%   threshold_proportional, which retains a FIXED PROPORTION of the strongest
%   edges. That operation depends only on the RANK ORDER of edge weights, so
%   any monotone transformation applied uniformly to every edge, additive or
%   multiplicative, leaves the binarized graph, and therefore the clustering
%   coefficient, exactly unchanged. A uniform inflation of connectivity
%   cannot by itself produce a clustering effect.
%
%   What proportional thresholding does NOT protect against is a NON-UNIFORM
%   perturbation, one that affects some channel pairs more than others, since
%   that does reorder edges. This function quantifies both cases on the real
%   connectivity matrices, giving the magnitude of each.
%
%   Four variants are compared, per subject, over the sparsity range:
%       1  baseline
%       2  every edge shifted by +delta            (uniform, additive)
%       3  every edge scaled by (1 + delta)        (uniform, multiplicative)
%       4  every edge perturbed by N(0, delta)     (non-uniform)
%
%   INPUT
%       dataDir     project 'data' folder
%       condition   'Lap' or 'Resting'
%       chromophore 'HbO', 'HbR' or 'HbT'
%       nSubUse     subjects to use (default all)
%       seed        rng seed (default 42)
%
%   OUTPUT struct I with the mean clustering coefficient under each variant
%   and the change relative to baseline.

if nargin < 2 || isempty(condition),   condition = 'Lap'; end
if nargin < 3 || isempty(chromophore), chromophore = 'HbO'; end
if nargin < 5 || isempty(seed),        seed = 42; end
rng(seed);

delta = 0.05;                       % comparable to the observed FC difference
threshold = 0.1:0.01:0.34;

[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];

f = fullfile(dataDir, sprintf('mid%sConn%s.mat', condition, chromophore));
assert(isfile(f), 'Not found: %s', f);
S = load(f, 'zMatFDR', 'keepRun');

nSub = numel(S.zMatFDR);
if nargin >= 4 && ~isempty(nSubUse), nSub = min(nSub, nSubUse); end

variants = {'baseline', 'uniform additive', 'uniform multiplicative', 'non-uniform'};
C  = nan(nSub, numel(variants));   % corrected, binarized clustering
CW = nan(nSub, numel(variants));   % as submitted, weighted

for iSub = 1:nSub
    z = S.zMatFDR{iSub};
    if isempty(z), continue; end
    Z = z(channelIdx, channelIdx);

    for v = 1:numel(variants)
        switch v
            case 1, W = Z;
            case 2, W = Z + delta;
            case 3, W = Z * (1 + delta);
            case 4
                E = delta * randn(size(Z));
                E = triu(E, 1); E = E + E';
                W = Z + E;
        end
        cc  = nan(1, numel(threshold));
        ccw = nan(1, numel(threshold));
        for iThr = 1:numel(threshold)
            Wt = threshold_proportional(W, threshold(iThr));
            Wb = weight_conversion(Wt, 'binarize');
            cc(iThr)  = mean(clustering_coef_bu(Wb), 'omitnan');   % corrected
            ccw(iThr) = mean(clustering_coef_bu(Wt), 'omitnan');   % as submitted
        end
        C(iSub, v)  = mean(cc(isfinite(cc)));
        CW(iSub, v) = mean(ccw(isfinite(ccw)));
    end
end

I.delta       = delta;
I.variants    = variants;
I.clustering   = C;
I.clusteringW  = CW;
I.meanClust    = mean(C,  1, 'omitnan');
I.meanClustW   = mean(CW, 1, 'omitnan');
I.deltaClust   = I.meanClust  - I.meanClust(1);
I.deltaClustW  = I.meanClustW - I.meanClustW(1);
I.pctChange    = 100 * I.deltaClust  / I.meanClust(1);
I.pctChangeW   = 100 * I.deltaClustW / I.meanClustW(1);

fprintf('\nClustering under simulated connectivity inflation (delta = %.2f)\n', delta);
fprintf('%-26s %12s %9s   %12s %9s\n', '', 'binarized', 'change', 'weighted', 'change');
for v = 1:numel(variants)
    fprintf('  %-24s %12.6f %+8.2f%%   %12.6f %+8.2f%%\n', variants{v}, ...
            I.meanClust(v), I.pctChange(v), I.meanClustW(v), I.pctChangeW(v));
end
fprintf(['\nBinarized, a uniform inflation leaves the metric unchanged, because\n' ...
         'proportional thresholding depends only on edge rank order; only the\n' ...
         'non-uniform perturbation moves it. Computed on the weighted graph, as\n' ...
         'in the submitted pipeline, the metric scales with the inflation itself.\n']);

end

% EOF
