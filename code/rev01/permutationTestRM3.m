function [p, observeddifference, Fobs] = permutationTestRM3(X, permutations, statistic)
% PERMUTATIONTESTRM3  Restricted permutation test for one within-subject factor.
%
%   [p, observeddifference, Fobs] = permutationTestRM3(X, permutations)
%   [p, observeddifference, Fobs] = permutationTestRM3(X, permutations, statistic)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   The submitted analysis used permutationTest3, which pools all 90
%   observations and permutes freely. That scheme assumes independent
%   samples, whereas the same 30 residents were measured at 0, 12 and 24 h.
%   Free permutation breaks the subject pairing and yields an anticonservative
%   null distribution. Here the three condition labels are permuted WITHIN
%   each subject, which is the exchangeability structure implied by the
%   repeated-measures design.
%
%   INPUT
%       X            [nSub x 3] matrix, one row per subject, one column per
%                    time point (0 h, 12 h, 24 h). NaNs are tolerated.
%       permutations scalar, number of permutations (e.g. 5000).
%       statistic    'maxmin' (default) reproduces the statistic used by
%                    permutationTest3, max(colmean) - min(colmean).
%                    'F' uses the repeated-measures F ratio, which is more
%                    powerful and is reported alongside either way.
%
%   OUTPUT
%       p                  permutation p-value (two-sided by construction)
%       observeddifference max(colmean) - min(colmean) in the observed data
%       Fobs               observed repeated-measures F ratio
%
%   Only rows with complete data across the three time points are used, so
%   that the within-subject permutation is well defined.

if nargin < 2 || isempty(permutations), permutations = 5000; end
if nargin < 3 || isempty(statistic),    statistic = 'maxmin'; end

% Complete cases only: within-subject permutation needs all three levels
keep = all(~isnan(X), 2);
X    = X(keep, :);
[nSub, nCond] = size(X);

if nCond ~= 3
    error('permutationTestRM3:badSize', 'X must have exactly three columns.');
end
if nSub < 3
    warning('permutationTestRM3:tooFewSubjects', ...
            'Only %d complete cases; returning NaN.', nSub);
    p = NaN; observeddifference = NaN; Fobs = NaN;
    return
end

    function [d, F] = stats(M)
        cm = mean(M, 1);
        d  = max(cm) - min(cm);
        % One-way repeated-measures ANOVA on M [nSub x nCond]
        gm      = mean(M(:));
        sm      = mean(M, 2);                       % subject means
        SS_cond = nSub * sum((cm - gm).^2);
        SS_subj = nCond * sum((sm - gm).^2);
        SS_tot  = sum((M(:) - gm).^2);
        SS_err  = SS_tot - SS_cond - SS_subj;
        df_cond = nCond - 1;
        df_err  = (nSub - 1) * (nCond - 1);
        F       = (SS_cond / df_cond) / (SS_err / df_err);
    end

[observeddifference, Fobs] = stats(X);

switch lower(statistic)
    case 'maxmin', obs = observeddifference;
    case 'f',      obs = Fobs;
    otherwise, error('permutationTestRM3:badStat', 'Unknown statistic ''%s''.', statistic);
end

permStat = nan(permutations, 1);
for iPerm = 1:permutations
    Xp = X;
    for iSub = 1:nSub
        Xp(iSub, :) = X(iSub, randperm(nCond));   % permute labels within subject
    end
    [d, F] = stats(Xp);
    if strcmpi(statistic, 'f'), permStat(iPerm) = F; else, permStat(iPerm) = d; end
end

% Both statistics are non-negative and large under the alternative,
% so the test is one-tailed on the statistic and two-sided on the effect.
% Add-one correction keeps p strictly positive (Phipson & Smyth, 2010).
p = (sum(permStat >= obs) + 1) / (permutations + 1);

end

% EOF
