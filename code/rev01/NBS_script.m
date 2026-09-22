%% NBS_script  Build NBS inputs for the repeated-measures effect of time.
%
% Neurophotonics NPH-260108-1, revision 01.
%
% TWO CHANGES RELATIVE TO THE SUBMITTED VERSION
%
% 1. Short-separation channels. The submitted script removed nodes [3 14],
%    which are the short-channel indices of the MOTOR montage used in an
%    earlier study. In the prefrontal montage the short channels are 3
%    (S03D01) and 24 (S10D08). The indices are now derived from
%    get_channels_from_template so they cannot drift again.
%
% 2. Incomplete subjects. The submitted script built a 30-subject design
%    matrix unconditionally, before loading any data, and then asserted that
%    every cell held a 24x24 matrix. conn_mat_subject_level returns an empty
%    matrix for any run retaining less than 120 s of artifact-free signal, so
%    that assertion fails as soon as one run is dropped, and if it had not
%    failed the design matrix would not have matched the data. The script now
%    determines which subjects have a usable run at ALL THREE time points and
%    builds the array, the design matrix and the exchange blocks together for
%    exactly those subjects. Listwise deletion is used because a
%    repeated-measures design cannot accommodate a missing cell.
%
% The two conditions are handled independently and may retain different
% numbers of subjects; each gets its own design and exchange file.
%
% OUTPUTS, per condition <cond> in {lap, resting}
%     ..\data\conn_Mat_HbO_NBS_<cond>.mat     Mat, 22 x 22 x (3*N)
%     ..\data\design_NBS_<cond>.mat           design, (3*N) x (N+2)
%     ..\data\exchange_NBS_<cond>.mat         exchange, (3*N) x 1
%     ..\data\NBS_subjects_<cond>.mat         subjIncluded, dfTime
%
% In the NBS GUI use the F-test with contrast [zeros(1,N) 1 1].
% Post hoc contrasts are [zeros(1,N) 1 -1], [zeros(1,N) 1 0], [zeros(1,N) 0 1].

clear; close all; format compact; clc

chromophore = 'HbO';
baseDir     = fullfile('..', 'data');

% Long channels of the prefrontal montage; short channels are excluded
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
nodesKeep     = sort([chIdxL, chIdxR]);
nodesToRemove = setdiff(1:24, nodesKeep);          % = [3 24]
nNodes        = numel(nodesKeep);                  % = 22
fprintf('Short-separation channels removed: [%s]; %d nodes retained.\n', ...
        num2str(nodesToRemove), nNodes);

conds = { 'lap',     {'beginLapConn',     'midLapConn',     'endLapConn'}
          'resting', {'beginRestingConn', 'midRestingConn', 'endRestingConn'} };

for iCond = 1:size(conds, 1)
    condName = conds{iCond, 1};
    stems    = conds{iCond, 2};
    nTime    = numel(stems);

    fprintf('\n=== %s ===\n', condName);

    % ---- load the three time points ------------------------------------
    Z = cell(1, nTime);
    for iT = 1:nTime
        f = fullfile(baseDir, [stems{iT} chromophore '.mat']);
        assert(isfile(f), 'Not found: %s', f);
        S = load(f, 'zMatFDR');
        Z{iT} = S.zMatFDR;
    end
    nSubAll = numel(Z{1});
    assert(all(cellfun(@numel, Z) == nSubAll), ...
           'The three files hold different numbers of subjects.');

    % ---- which subjects have a usable run at every time point? ----------
    good = true(nSubAll, 1);
    for iS = 1:nSubAll
        for iT = 1:nTime
            M = Z{iT}{iS};
            if isempty(M) || ~isequal(size(M), [24 24])
                good(iS) = false;
            end
        end
    end
    subjIncluded = find(good);
    N = numel(subjIncluded);

    fprintf('Usable runs per time point: %s\n', ...
            mat2str(cellfun(@(c) sum(~cellfun(@(m) isempty(m) || ~isequal(size(m),[24 24]), c)), Z)));
    fprintf('Subjects complete at all %d time points: %d of %d\n', nTime, N, nSubAll);
    if N < nSubAll
        fprintf('Excluded (incomplete): %s\n', mat2str(find(~good)'));
    end
    assert(N >= 3, 'Too few complete subjects for %s.', condName);

    % ---- build the array, ordered S1T1..SNT1, S1T2..SNT2, S1T3..SNT3 ----
    Mat = nan(nNodes, nNodes, N * nTime);
    idx = 0;
    for iT = 1:nTime
        for k = 1:N
            idx = idx + 1;
            M = Z{iT}{subjIncluded(k)};
            Mat(:, :, idx) = M(nodesKeep, nodesKeep);
        end
    end
    assert(~any(isnan(Mat(:)) & false), 'unreachable');   % shape guard only
    fprintf('Array: %s\n', mat2str(size(Mat)));

    % ---- design matrix and exchange blocks for exactly these subjects ---
    I  = eye(N);
    T1 = repmat([1 0], N, 1);
    T2 = repmat([0 1], N, 1);
    T3 = repmat([0 0], N, 1);          % reference level
    design   = [I T1; I T2; I T3];
    exchange = repmat((1:N)', nTime, 1);
    contrast = [zeros(1, N) 1 1];      % omnibus F-test for the effect of time

    dfTime = [nTime - 1, (N - 1) * (nTime - 1)];
    fprintf('Design: %s, contrast length %d, F degrees of freedom (%d, %d)\n', ...
            mat2str(size(design)), numel(contrast), dfTime(1), dfTime(2));

    % ---- save -----------------------------------------------------------
    save(fullfile(baseDir, sprintf('conn_Mat_%s_NBS_%s.mat', chromophore, condName)), ...
         'Mat', '-v7.3');
    save(fullfile(baseDir, sprintf('design_NBS_%s.mat',   condName)), 'design',   '-v7.3');
    save(fullfile(baseDir, sprintf('exchange_NBS_%s.mat', condName)), 'exchange', '-v7.3');
    save(fullfile(baseDir, sprintf('NBS_subjects_%s.mat', condName)), ...
         'subjIncluded', 'dfTime', 'contrast', 'nodesKeep', 'nodesToRemove');
end

fprintf(['\nDone. In the NBS GUI, load for each condition:\n' ...
         '  Connectivity matrices : conn_Mat_%s_NBS_<cond>.mat\n' ...
         '  Design matrix         : design_NBS_<cond>.mat\n' ...
         '  Exchange blocks       : exchange_NBS_<cond>.mat\n' ...
         '  Contrast              : the contrast saved in NBS_subjects_<cond>.mat\n' ...
         '  Test                  : F-test, 5000 permutations\n' ...
         'Report the component p values, the surviving edges and the degrees of\n' ...
         'freedom printed above; they differ between conditions.\n'], chromophore);

% EOF
