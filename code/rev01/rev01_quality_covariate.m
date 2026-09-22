function Q = rev01_quality_covariate(M, SCIsubj, PSPsubj, pctMotion, condition, nPerms, seed)
% REV01_QUALITY_COVARIATE  Does the time effect survive adjustment for data quality?
%
%   Q = rev01_quality_covariate(M, SCIsubj, PSPsubj, pctMotion, condition, nPerms, seed)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Written for Associate Editor major concern 1d. The mid-shift clustering
%   effect coincides with the session that scores worst on every quality
%   metric and that also shows the highest mean connectivity. Poor scalp
%   coupling and residual motion can both inflate correlation, and hence
%   density-dependent clustering, so time and data quality are entangled.
%   This function quantifies that entanglement and tests whether the effect
%   of time survives once quality is accounted for.
%
%   Three complementary analyses are run:
%     A. Descriptive. Mean SCI, PSP, motion percentage and mean connectivity
%        at each time point, with the repeated-measures test for each.
%     B. Linear mixed-effects models with a random intercept per subject:
%           AUCclust ~ 1 + Time                      (unadjusted)
%           AUCclust ~ 1 + Time + SCI + pctMotion    (adjusted)
%        The effect of Time is tested in both by Type III ANOVA, and the two
%        models are compared by likelihood ratio.
%     C. Permutation. The repeated-measures permutation test of the
%        submitted framework, applied to the residuals of
%        AUCclust ~ SCI + pctMotion, so that the resampling test is run on
%        quality-adjusted data.
%
%   INPUT
%       M          struct returned by rev01_network_measures
%       SCIsubj    [nSub x 6] per-scan mean SCI      from qc_summary_export
%       PSPsubj    [nSub x 6] per-scan mean PSP      from qc_summary_export
%       pctMotion  [nSub x 6] percentage of time flagged as motion artifact,
%                             from motion_artifacts_summary.mat
%       condition  'Lap' or 'Resting'; selects columns 1:3 or 4:6
%       nPerms     permutations (default 5000)
%       seed       rng seed (default 42)
%
%   Column order of SCIsubj, PSPsubj and pctMotion must be
%   {lap01 lap02 lap03 resting01 resting02 resting03}.

if nargin < 6 || isempty(nPerms), nPerms = 5000; end
if nargin < 7 || isempty(seed),   seed   = 42;   end
rng(seed);

switch lower(condition)
    case 'lap',     cols = 1:3;
    case 'resting', cols = 4:6;
    otherwise, error('rev01_quality_covariate:badCondition', ...
                     'condition must be ''Lap'' or ''Resting''.');
end

SCI    = SCIsubj(:, cols);
PSP    = PSPsubj(:, cols);
MOT    = pctMotion(:, cols);
CLUST  = M.AUCclust;
LAMBDA = M.AUClambda;
SWIa   = M.AUCswi;
FC     = M.meanFC;

nSub = size(CLUST, 1);
if size(SCI, 1) ~= nSub
    error('rev01_quality_covariate:sizeMismatch', ...
          'Quality matrices have %d rows, network metrics have %d.', size(SCI,1), nSub);
end

Q.condition = condition;

%% A. Descriptives and repeated-measures tests on the quality metrics
Q.meanSCI = mean(SCI, 1, 'omitnan');
Q.meanPSP = mean(PSP, 1, 'omitnan');
Q.meanMOT = mean(MOT, 1, 'omitnan');
Q.meanFC  = mean(FC,  1, 'omitnan');
Q.sdSCI   = std(SCI, 0, 1, 'omitnan');
Q.sdPSP   = std(PSP, 0, 1, 'omitnan');
Q.sdMOT   = std(MOT, 0, 1, 'omitnan');
Q.sdFC    = std(FC,  0, 1, 'omitnan');

Q.pSCI = permutationTestRM3(SCI, nPerms);
Q.pPSP = permutationTestRM3(PSP, nPerms);
Q.pMOT = permutationTestRM3(MOT, nPerms);
Q.pFC  = permutationTestRM3(FC,  nPerms);

%% Association of quality with connectivity and with clustering, pooled
ok = ~isnan(SCI(:)) & ~isnan(MOT(:)) & ~isnan(FC(:)) & ~isnan(CLUST(:));
[Q.rSCI_FC,    Q.pSCI_FC]    = corr(SCI(ok),   FC(ok));
[Q.rMOT_FC,    Q.pMOT_FC]    = corr(MOT(ok),   FC(ok));
[Q.rSCI_CLUST, Q.pSCI_CLUST] = corr(SCI(ok),   CLUST(ok));
[Q.rMOT_CLUST, Q.pMOT_CLUST] = corr(MOT(ok),   CLUST(ok));
[Q.rFC_CLUST,  Q.pFC_CLUST]  = corr(FC(ok),    CLUST(ok));

%% B. Linear mixed-effects models
subject = repmat((1:nSub)', 3, 1);
timeNum = reshape(repmat(1:3, nSub, 1), [], 1);
timeLab = categorical(timeNum, 1:3, {'h00', 'h12', 'h24'});

T = table(categorical(subject), timeLab, CLUST(:), LAMBDA(:), SWIa(:), ...
          FC(:), SCI(:), PSP(:), MOT(:), ...
          'VariableNames', {'Subject', 'Time', 'AUCclust', 'AUClambda', ...
                            'AUCswi', 'meanFC', 'SCI', 'PSP', 'pctMotion'});
T = rmmissing(T);
Q.table = T;
Q.nObsUsed = height(T);

Q.lmeOK = true;
try
    Q.lmeUnadj = fitlme(T, 'AUCclust ~ 1 + Time + (1|Subject)');
    Q.lmeAdj   = fitlme(T, 'AUCclust ~ 1 + Time + SCI + pctMotion + (1|Subject)');
    Q.anovaUnadj = anova(Q.lmeUnadj, 'DFMethod', 'satterthwaite');
    Q.anovaAdj   = anova(Q.lmeAdj,   'DFMethod', 'satterthwaite');
catch ME
    warning('rev01_quality_covariate:lmeFailed', ...
            ['fitlme failed (%s). Falling back to Satterthwaite-free ' ...
             'default degrees of freedom.'], ME.message);
    Q.lmeUnadj = fitlme(T, 'AUCclust ~ 1 + Time + (1|Subject)');
    Q.lmeAdj   = fitlme(T, 'AUCclust ~ 1 + Time + SCI + pctMotion + (1|Subject)');
    Q.anovaUnadj = anova(Q.lmeUnadj);
    Q.anovaAdj   = anova(Q.lmeAdj);
end

% Extract the Time row of each Type III ANOVA table
Q.pTimeUnadj = Q.anovaUnadj.pValue(strcmp(Q.anovaUnadj.Term, 'Time'));
Q.FTimeUnadj = Q.anovaUnadj.FStat( strcmp(Q.anovaUnadj.Term, 'Time'));
Q.pTimeAdj   = Q.anovaAdj.pValue(  strcmp(Q.anovaAdj.Term,   'Time'));
Q.FTimeAdj   = Q.anovaAdj.FStat(   strcmp(Q.anovaAdj.Term,   'Time'));
Q.dfTimeAdj  = [Q.anovaAdj.DF1(strcmp(Q.anovaAdj.Term,'Time')), ...
                Q.anovaAdj.DF2(strcmp(Q.anovaAdj.Term,'Time'))];

% Covariate coefficients in the adjusted model
cf = Q.lmeAdj.Coefficients;
Q.betaSCI   = cf.Estimate(strcmp(cf.Name, 'SCI'));
Q.pBetaSCI  = cf.pValue(  strcmp(cf.Name, 'SCI'));
Q.betaMOT   = cf.Estimate(strcmp(cf.Name, 'pctMotion'));
Q.pBetaMOT  = cf.pValue(  strcmp(cf.Name, 'pctMotion'));

% Likelihood-ratio comparison of the two models
Q.lrt = compare(Q.lmeUnadj, Q.lmeAdj, 'CheckNesting', true);

%% C. Permutation test on quality-adjusted residuals
% Residualize AUCclust on SCI and pctMotion, ignoring time, then apply the
% same repeated-measures permutation test used for the unadjusted metric.
X = [ones(height(T),1), T.SCI, T.pctMotion];
b = X \ T.AUCclust;
res = T.AUCclust - X * b;

R = nan(nSub, 3);
subjNum = double(string(T.Subject));
timeIdx = double(T.Time);
for k = 1:height(T)
    R(subjNum(k), timeIdx(k)) = res(k);
end
Q.residClust = R;
[Q.pResidClust, Q.dResidClust] = permutationTestRM3(R, nPerms);

pairs = {[1 2], [1 3], [2 3]};
names = {'12', '13', '23'};
for k = 1:3
    a = pairs{k}(1); b2 = pairs{k}(2);
    Q.(['pResidClust_' names{k}]) = ...
        permutationTestRMpair(R(:, a), R(:, b2), nPerms);
end

end

% EOF
