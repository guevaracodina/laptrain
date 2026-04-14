function [p, observeddifference, effectsize] = permutationTest3(sample1, sample2, sample3, permutations, varargin)
% [p, observeddifference, effectsize] = permutationTest3(sample1, sample2, sample3, permutations [, varargin])
%
% Permutation test for three independent samples.
% Tests for a difference in means between three samples.
%
% Inputs:
%   sample1, sample2, sample3 - vectors of measurements
%   permutations - number of permutations
%
% Optional (name-value pairs):
%   sidedness - 'both' (default), 'smaller', 'larger'
%   exact - 0 (default) or 1; whether to perform an exact test
%   plotresult - 0 (default) or 1; plot distribution of permuted differences
%   showprogress - 0 (default); if >0, update progress bar every N iterations
%
% Outputs:
%   p - p-value
%   observeddifference - observed statistic (max(mean) - min(mean))
%   effectsize - generalized eta squared (Cohen 1988)
%
% Example:
%   permutationTest3(randn(1,30), randn(1,30)+1, randn(1,30)+2, 10000, 'plotresult', 1)

% Parse inputs
p = inputParser;
addRequired(p, 'sample1', @isnumeric);
addRequired(p, 'sample2', @isnumeric);
addRequired(p, 'sample3', @isnumeric);
addRequired(p, 'permutations', @isnumeric);
addParameter(p, 'sidedness', 'both', @(x) any(validatestring(x,{'both', 'smaller', 'larger'})));
addParameter(p, 'exact' , 0, @isnumeric);
addParameter(p, 'plotresult', 0, @isnumeric);
addParameter(p, 'showprogress', 0, @isnumeric);
parse(p, sample1, sample2, sample3, permutations, varargin{:});

sample1 = p.Results.sample1;
sample2 = p.Results.sample2;
sample3 = p.Results.sample3;
permutations = p.Results.permutations;
sidedness = p.Results.sidedness;
exact = p.Results.exact;
plotresult = p.Results.plotresult;
showprogress = p.Results.showprogress;

% Ensure row vectors
if iscolumn(sample1), sample1 = sample1'; end
if iscolumn(sample2), sample2 = sample2'; end
if iscolumn(sample3), sample3 = sample3'; end

allobservations = [sample1, sample2, sample3];
n1 = numel(sample1);
n2 = numel(sample2);
n3 = numel(sample3);

% Test statistic: max mean - min mean
means = [nanmean(sample1), nanmean(sample2), nanmean(sample3)];
observeddifference = max(means) - min(means);

% Generalized eta squared (see: Cohen, 1988)
groupmeans = means;
grandmean = nanmean(allobservations);
SStotal = nansum((allobservations - grandmean).^2);
SSbetween = n1*(groupmeans(1)-grandmean)^2 + n2*(groupmeans(2)-grandmean)^2 + n3*(groupmeans(3)-grandmean)^2;
effectsize = SSbetween/SStotal; % Generalized eta squared

if showprogress, w = waitbar(0, 'Preparing test...', 'Name', 'permutationTest3'); end

if exact
    error('Exact test not implemented for three samples due to combinatorial explosion.');
end

randomdifferences = zeros(1, permutations);
if showprogress, waitbar(0, w, sprintf('Permutation 1 of %d', permutations), 'Name', 'permutationTest3'); end
for n = 1:permutations
    if showprogress && mod(n,showprogress) == 0
        waitbar(n/permutations, w, sprintf('Permutation %d of %d', n, permutations));
    end
    perm = randperm(length(allobservations));
    randSample1 = allobservations(perm(1:n1));
    randSample2 = allobservations(perm(n1+1:n1+n2));
    randSample3 = allobservations(perm(n1+n2+1:end));
    randommeans = [nanmean(randSample1), nanmean(randSample2), nanmean(randSample3)];
    randomdifferences(n) = max(randommeans) - min(randommeans);
end
if showprogress, delete(w); end

% Compute p-value
switch sidedness
    case 'both'
        p = (sum(abs(randomdifferences) >= abs(observeddifference))+1) / (permutations+1);
    case 'smaller'
        p = (sum(randomdifferences < observeddifference)+1) / (permutations+1);
    case 'larger'
        p = (sum(randomdifferences > observeddifference)+1) / (permutations+1);
end

% Plot result
if plotresult
    figure;
    histogram(randomdifferences, 20);
    hold on;
    xlabel('Random differences');
    ylabel('Count')
    od = plot(observeddifference, 0, '*r', 'DisplayName', sprintf('Observed difference.\nEffect size: %.2f,\np = %f', effectsize, p));
    legend(od);
end

end
