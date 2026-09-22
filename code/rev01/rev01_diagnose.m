function D = rev01_diagnose(dataDir)
% REV01_DIAGNOSE  Two questions that must be settled before the revision proceeds.
%
%   D = rev01_diagnose(dataDir)
%
%   Neurophotonics NPH-260108-1, revision 01.
%
%   QUESTION 1. Which data column is which channel?
%   get_channels_from_template('prefrontal') lists the montage ordered by
%   DETECTOR, putting the short channels at indices 3 and 24. The corrected
%   optodes_MNI_coordinates.csv lists it ordered by SOURCE, which puts S03D01
%   at row 5 and S10D08 at row 24. Both name the same two short channels, but
%   they disagree about the INDEX those channels occupy. Everything downstream
%   depends on which ordering the Homer3 output actually uses. This function
%   reads it directly from output.dc.measurementList and says so.
%
%   QUESTION 2. Is the clustering coefficient being computed on a binary graph?
%   clustering_coef_bu expects a BINARY undirected matrix; it computes
%   sum(S(:))/(k^2-k) over the neighbour submatrix. Passing it the weighted
%   output of threshold_proportional makes that sum a sum of WEIGHTS, so the
%   result scales linearly with connectivity strength instead of measuring
%   topology. This function scales every edge by a known factor and reports
%   whether the metric scales with it.
%
%   INPUT
%       dataDir  project 'data' folder. Default '..\data'.

if nargin < 1 || isempty(dataDir), dataDir = fullfile('..', 'data'); end

D = struct();

%% ---------------------------------------------------------------- Q1 ----
fprintf('\n%s\n QUESTION 1: actual channel order in the Homer3 output\n%s\n', ...
        repmat('=',1,72), repmat('=',1,72));

f = fullfile(dataDir, 'lap01', 'derivatives', 'homer', 'sub-01', 'nirs', ...
             'sub-01_task-lap01_nirs.mat');
if ~isfile(f)
    error('rev01_diagnose:noFile', ['Not found: %s\nPoint dataDir at the ' ...
          'folder holding the extracted derivatives.'], f);
end
S = load(f, 'output');
ml = S.output.dc.measurementList;

srcList = []; detList = [];
for k = 1:numel(ml)
    if contains(ml(k).dataTypeLabel, 'HbO')
        srcList(end+1) = ml(k).sourceIndex;   %#ok<AGROW>
        detList(end+1) = ml(k).detectorIndex; %#ok<AGROW>
    end
end
D.sourceIndex = srcList(:);
D.detectorIndex = detList(:);
D.nChannels = numel(srcList);

fprintf('\nHbO channels found: %d\n\n  idx   S-D pair\n', D.nChannels);
for k = 1:D.nChannels
    fprintf('  %3d   S%02dD%02d\n', k, srcList(k), detList(k));
end

isSSC = (srcList == 3 & detList == 1) | (srcList == 10 & detList == 8);
D.sscIdxFromData = find(isSSC);
fprintf('\n  S03D01 and S10D08 sit at data indices: [%s]\n', num2str(D.sscIdxFromData));

detSorted = all(diff(detList) >= 0);
srcSorted = all(diff(srcList) >= 0);
fprintf('  Ordering is detector-major: %d    source-major: %d\n', detSorted, srcSorted);

if isequal(D.sscIdxFromData(:)', [3 24])
    fprintf('\n  >> The template is correct. SSC = [3 24]. The re-run stands.\n');
elseif isequal(D.sscIdxFromData(:)', [5 24])
    fprintf(['\n  >> The corrected CSV ordering is correct. SSC = [5 24].\n' ...
             '     conn_regress_global and NBS_script must be re-run with [5 24],\n' ...
             '     and get_channels_from_template must be rewritten.\n']);
else
    fprintf('\n  >> Neither expected pattern. Indices are [%s].\n', num2str(D.sscIdxFromData));
end

%% ---------------------------------------------------------------- Q2 ----
fprintf('\n%s\n QUESTION 2: does the clustering coefficient scale with weights?\n%s\n', ...
        repmat('=',1,72), repmat('=',1,72));

fc = fullfile(dataDir, 'midLapConnHbO.mat');
assert(isfile(fc), 'Not found: %s', fc);
C = load(fc, 'zMatFDR');

[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
channelIdx = [chIdxL chIdxR];
threshold = 0.1:0.01:0.34;
scales = [1.00 1.05 1.25 2.00];

nUse = min(10, numel(C.zMatFDR));
weighted = nan(nUse, numel(scales));
binary   = nan(nUse, numel(scales));

for iSub = 1:nUse
    z = C.zMatFDR{iSub};
    if isempty(z), continue; end
    Z = z(channelIdx, channelIdx);
    for s = 1:numel(scales)
        vw = nan(1, numel(threshold));
        vb = nan(1, numel(threshold));
        for i = 1:numel(threshold)
            Wt = threshold_proportional(Z * scales(s), threshold(i));
            vw(i) = mean(clustering_coef_bu(Wt), 'omitnan');                    % as used now
            vb(i) = mean(clustering_coef_bu(weight_conversion(Wt,'binarize')), 'omitnan');
        end
        weighted(iSub, s) = mean(vw(isfinite(vw)));
        binary(iSub, s)   = mean(vb(isfinite(vb)));
    end
end

D.scales = scales;
D.clustWeighted = mean(weighted, 1, 'omitnan');
D.clustBinary   = mean(binary,   1, 'omitnan');
D.ratioWeighted = D.clustWeighted / D.clustWeighted(1);
D.ratioBinary   = D.clustBinary   / D.clustBinary(1);

fprintf('\n  scale   C as currently computed   ratio      C if binarized   ratio\n');
for s = 1:numel(scales)
    fprintf('  x%.2f   %18.6f   %.4f   %14.6f   %.4f\n', scales(s), ...
            D.clustWeighted(s), D.ratioWeighted(s), D.clustBinary(s), D.ratioBinary(s));
end

if max(abs(D.ratioWeighted - scales)) < 0.01
    fprintf(['\n  >> The metric scales EXACTLY with the edge weights. It is not a\n' ...
             '     binary clustering coefficient; it is a weight sum, so it tracks\n' ...
             '     mean connectivity by construction.\n']);
else
    fprintf('\n  >> The metric does not scale exactly with weights.\n');
end

end

% EOF
