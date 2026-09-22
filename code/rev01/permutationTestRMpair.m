function [p, observeddifference, dz] = permutationTestRMpair(x, y, permutations)
% PERMUTATIONTESTRMPAIR  Exact-style sign-flip permutation test for paired data.
%
%   [p, observeddifference, dz] = permutationTestRMpair(x, y, permutations)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Post-hoc counterpart of permutationTestRM3. The submitted analysis used
%   permutationTest, which permutes group membership as though the two
%   samples were independent. For measurements repeated in the same subject
%   the exchangeable quantity is the SIGN of the within-subject difference,
%   so the null distribution is built by randomly flipping those signs.
%
%   INPUT
%       x, y         vectors of equal length, paired observations
%       permutations scalar, number of sign-flip permutations
%
%   OUTPUT
%       p                  permutation p-value, two-sided
%       observeddifference mean(x - y)
%       dz                 Cohen's d for paired samples, mean(d)/std(d)

if nargin < 3 || isempty(permutations), permutations = 5000; end

x = x(:); y = y(:);
if numel(x) ~= numel(y)
    error('permutationTestRMpair:badSize', 'x and y must be the same length.');
end

d = x - y;
d = d(~isnan(d));
n = numel(d);
if n < 3
    warning('permutationTestRMpair:tooFewPairs', ...
            'Only %d complete pairs; returning NaN.', n);
    p = NaN; observeddifference = NaN; dz = NaN;
    return
end

observeddifference = mean(d);
dz = mean(d) / std(d);

permStat = nan(permutations, 1);
for iPerm = 1:permutations
    signs = 2 * (rand(n, 1) > 0.5) - 1;       % random +1/-1
    permStat(iPerm) = mean(signs .* d);
end

p = (sum(abs(permStat) >= abs(observeddifference)) + 1) / (permutations + 1);

end

% EOF
