function rev01_write_report(R, outFile)
% REV01_WRITE_REPORT  Human-readable summary of the revision-01 analyses.
%
%   rev01_write_report(R, outFile)
%
%   R is the struct assembled by rev01_run_all. Everything needed to write
%   the response letter is printed here in one place.

REPORT_VERSION = '2026-09-18b';
fprintf('rev01_write_report version %s\n', REPORT_VERSION);

fid = fopen(outFile, 'w');
assert(fid > 0, 'Could not open %s for writing.', outFile);
c = onCleanup(@() fclose(fid));

p = @(varargin) fprintf(fid, varargin{:});

p('Neurophotonics NPH-260108-1, revision 01\n');
p('Report writer version %s\n', REPORT_VERSION);
p('Analysis report generated %s\n', R.generated);
p('Permutations: %d, rng seed: %d\n', R.nPerms, R.seed);
p('Short-separation channels: [%s]; %d long channels retained\n', ...
  num2str(R.sscIndices), numel(R.longChannels));
p('%s\n', repmat('=', 1, 78));

%% Network metrics
try

    p('\n1. NETWORK METRICS, AUC over sparsity 0.10 to 0.34\n');
    p('   Repeated-measures permutation test, labels permuted within subject.\n');
    fn = fieldnames(R.NM);
    for k = 1:numel(fn)
        M = R.NM.(fn{k});
        p('\n   --- %s, %s ---\n', M.condition, M.chromophore);
        p('   n complete cases: %d of %d\n', sum(all(~isnan(M.AUCclust), 2)), size(M.AUCclust,1));
        mets = {'clust', 'clustW', 'lambda', 'lambdaB', 'swi'};
        labs = {'Clustering coeff (binary)', 'Clustering coeff (as submitted)', ...
                'Char. path length (weighted)', 'Char. path length (binary)', ...
                'Small-world index'};
        for m = 1:numel(mets)
            need = {['AUC' mets{m}], ['pAUC' mets{m}], ['FAUC' mets{m}], ...
                    ['pAUC' mets{m} '_12'], ['pAUC' mets{m} '_13'], ['pAUC' mets{m} '_23']};
            if ~all(isfield(M, need))
                p('   %-28s  [not available in this results file]\n', labs{m});
                continue
            end
            A  = M.(['AUC' mets{m}]);
            pv = M.(['pAUC' mets{m}]);
            F  = M.(['FAUC' mets{m}]);
            p('   %-28s  AUC mean 0 h / 12 h / 24 h = %.4f / %.4f / %.4f\n', ...
              labs{m}, mean(A(:,1),'omitnan'), mean(A(:,2),'omitnan'), mean(A(:,3),'omitnan'));
            p('   %-28s  global p = %.4f  (F = %.3f)\n', '', pv, F);
            p('   %-28s  post hoc 0-12 p = %.4f, 0-24 p = %.4f, 12-24 p = %.4f\n', '', ...
              M.(['pAUC' mets{m} '_12']), M.(['pAUC' mets{m} '_13']), M.(['pAUC' mets{m} '_23']));
        end
        p('   Mean FC (Fisher z) 0 h / 12 h / 24 h = %.4f / %.4f / %.4f\n', ...
          mean(M.meanFC(:,1),'omitnan'), mean(M.meanFC(:,2),'omitnan'), mean(M.meanFC(:,3),'omitnan'));
        p('   Corr(mean FC, binary clustering)       = %.3f / %.3f / %.3f\n', M.rFCclust);

        if isfield(M, 'rFCclustW')

            p('   Corr(mean FC, clustering as submitted) = %.3f / %.3f / %.3f\n', M.rFCclustW);

        end
        nSig = sum(M.pThrClust <= 0.05);
        p('   Thresholds with a significant BINARY clustering effect: %d of %d\n', ...
          nSig, numel(M.pThrClust));

        if isfield(M, 'pThrClustW')

            p('   Thresholds significant with the submitted metric:       %d of %d\n', ...
              sum(M.pThrClustW <= 0.05), numel(M.pThrClustW));

        end
        if isfield(M, 'covClust')
            p('   Threshold coverage (metric defined) clust / lambda / SWI = %.2f / %.2f / %.2f\n', ...
              mean(M.covClust(:)), mean(M.covLambda(:)), mean(M.covSWI(:)));
        end
    end
catch ME_sec
    p('\n   [section could not be written: %s]\n', ME_sec.message);
    warning('rev01_write_report:section', ...
            'Section 1 network metrics failed: %s', ME_sec.message);
end

%% Quality covariates
try

    p('\n%s\n', repmat('=', 1, 78));
    p('\n2. DATA QUALITY AS A CONFOUND (Associate Editor major concern 1d)\n');
    fn = fieldnames(R.QA);
    for k = 1:numel(fn)
        Q = R.QA.(fn{k});
        p('\n   --- %s ---\n', Q.condition);
        p('   Observations entering the models: %d\n', Q.nObsUsed);
        p('   Mean SCI        0 h / 12 h / 24 h = %.3f / %.3f / %.3f   (p = %.4f)\n', Q.meanSCI, Q.pSCI);
        p('   Mean PSP        0 h / 12 h / 24 h = %.4f / %.4f / %.4f  (p = %.4f)\n', Q.meanPSP, Q.pPSP);
        p('   Motion time %%   0 h / 12 h / 24 h = %.2f / %.2f / %.2f  (p = %.4f)\n', Q.meanMOT, Q.pMOT);
        p('   Mean FC (z)     0 h / 12 h / 24 h = %.4f / %.4f / %.4f  (p = %.4f)\n', Q.meanFC, Q.pFC);
        p('\n   Quality-connectivity associations, pooled over time points:\n');
        p('     SCI    vs mean FC        r = %+.3f (p = %.4g)\n', Q.rSCI_FC, Q.pSCI_FC);
        p('     motion vs mean FC        r = %+.3f (p = %.4g)\n', Q.rMOT_FC, Q.pMOT_FC);
        p('     SCI    vs AUC clustering r = %+.3f (p = %.4g)\n', Q.rSCI_CLUST, Q.pSCI_CLUST);
        p('     motion vs AUC clustering r = %+.3f (p = %.4g)\n', Q.rMOT_CLUST, Q.pMOT_CLUST);
        p('     meanFC vs AUC clustering r = %+.3f (p = %.4g)\n', Q.rFC_CLUST, Q.pFC_CLUST);
        p('\n   Linear mixed-effects models, random intercept per subject:\n');
        p('     AUCclust ~ Time                    : F = %.3f, p = %.4f\n', Q.FTimeUnadj, Q.pTimeUnadj);
        p('     AUCclust ~ Time + SCI + pctMotion  : F = %.3f, p = %.4f, df = (%d, %.1f)\n', ...
          Q.FTimeAdj, Q.pTimeAdj, Q.dfTimeAdj(1), Q.dfTimeAdj(2));
        p('     beta SCI       = %+.4f (p = %.4g)\n', Q.betaSCI, Q.pBetaSCI);
        p('     beta pctMotion = %+.4f (p = %.4g)\n', Q.betaMOT, Q.pBetaMOT);
        p('\n   Permutation test on quality-adjusted residuals:\n');
        p('     global p = %.4f; post hoc 0-12 p = %.4f, 0-24 p = %.4f, 12-24 p = %.4f\n', ...
          Q.pResidClust, Q.pResidClust_12, Q.pResidClust_13, Q.pResidClust_23);
        p('\n   Type III ANOVA, adjusted model:\n');
        printCell(fid, Q.anovaAdjTable, '     ');
        p('\n   Likelihood-ratio comparison:\n');
        printCell(fid, Q.lrtTable, '     ');
    end
catch ME_sec
    p('\n   [section could not be written: %s]\n', ME_sec.message);
    warning('rev01_write_report:section', ...
            'Section 2 quality covariates failed: %s', ME_sec.message);
end

%% Threshold invariance
try

    if isfield(R, 'IV')
        p('\n%s\n', repmat('=', 1, 78));
        p('\n2b. SENSITIVITY OF CLUSTERING TO GLOBAL CONNECTIVITY INFLATION\n');
        p('    Proportional thresholding fixes edge density and depends only on\n');
        p('    edge rank order, so a uniform inflation cannot change the metric.\n');
        fn = fieldnames(R.IV);
        for k = 1:numel(fn)
            I = R.IV.(fn{k});
            p('\n    --- %s, delta = %.2f ---\n', fn{k}, I.delta);
            p('      %-24s %12s %9s   %12s %9s\n', '', 'binarized', 'change', 'weighted', 'change');

            for v = 1:numel(I.variants)

                if isfield(I, 'meanClustW')

                    p('      %-24s %12.6f %+8.2f%%   %12.6f %+8.2f%%\n', I.variants{v}, ...
                      I.meanClust(v), I.pctChange(v), I.meanClustW(v), I.pctChangeW(v));

                else

                    p('      %-24s %12.6f %+8.2f%%\n', I.variants{v}, ...
                      I.meanClust(v), I.pctChange(v));

                end
            end
        end
    end
catch ME_sec
    p('\n   [section could not be written: %s]\n', ME_sec.message);
    warning('rev01_write_report:section', ...
            'Section 2b threshold invariance failed: %s', ME_sec.message);
end

%% Pruning
try

    p('\n%s\n', repmat('=', 1, 78));
    p('\n3. CHANNEL PRUNING AND RUN EXCLUSION (Reviewer 1)\n');
    PR = R.PR;
    p('   Minimum retained recording time applied: %d s\n\n', PR.minRecordingTime);
    p('   %-12s %-22s %-10s %-16s\n', 'condition', 'channels pruned', 'max', 'runs below min');
    for iCond = 1:numel(PR.condNames)
        p('   %-12s %.1f (%.1f)%-12s %-10d %-16d\n', PR.condNames{iCond}, ...
          PR.meanChanPruned(iCond), PR.sdChanPruned(iCond), '', ...
          PR.maxChanPruned(iCond), PR.nRunsExcluded(iCond));
    end
    p('\n   Mean retained recording time per condition (s): %s\n', ...
      num2str(PR.meanRetainedSec, '%.1f  '));
    if isempty(PR.subjectsFullyExcluded)
        p('   No subject was excluded in every condition.\n');
    else
        p('   Subjects excluded in every condition: %s\n', num2str(PR.subjectsFullyExcluded));
    end
catch ME_sec
    p('\n   [section could not be written: %s]\n', ME_sec.message);
    warning('rev01_write_report:section', ...
            'Section 3 pruning failed: %s', ME_sec.message);
end

%% Participants and behaviour
try

    p('\n%s\n', repmat('=', 1, 78));
    p('\n4. PARTICIPANTS AND BEHAVIOUR (Associate Editor minor concerns 1 and 3)\n');
    BD = R.BD;
    if isfield(BD, 'n'), p('   N = %d\n', BD.n); end
    if isfield(BD, 'ageMean')
        p('   Age: %.1f (%.1f) years, range %d to %d\n', ...
          BD.ageMean, BD.ageSD, BD.ageMin, BD.ageMax);
    end
    if isfield(BD, 'sexCategories')
        p('   Sex: %s\n', catLine(BD.sexCategories, BD.sexCounts));
    end
    if isfield(BD, 'handCategories')
        p('   Handedness: %s\n', catLine(BD.handCategories, BD.handCounts));
    end
    if isfield(BD, 'missingDemographics') && ~isempty(BD.missingDemographics)
        p('   NOT PRESENT in participants.tsv: %s\n', strjoin(BD.missingDemographics, ', '));
    end

    p('\n   Behaviour\n');
    if isfield(BD, 'PSCMean')
        p('     PSC          0 h / 12 h / 24 h = %.2f (%.2f) / %.2f (%.2f) / %.2f (%.2f)\n', ...
          BD.PSCMean(1), BD.PSCSD(1), BD.PSCMean(2), BD.PSCSD(2), BD.PSCMean(3), BD.PSCSD(3));
        p('     PSC          global p = %.4f; 0-12 p = %.4f, 0-24 p = %.4f, 12-24 p = %.4f\n', ...
          BD.PSCP, BD.PSCP12, BD.PSCP13, BD.PSCP23);
    end
    if isfield(BD, 'GRSConstantAcrossTime') && BD.GRSConstantAcrossTime
        p('     GRS          IDENTICAL at all three time points for every subject.\n');
        p('                  Single between-subject score: %.2f (%.2f), range %d to %d.\n', ...
          BD.GRSMean, BD.GRSSD, BD.GRSMin, BD.GRSMax);
        p('                  It cannot be reported as a repeated measure.\n');
    elseif isfield(BD, 'GRSMean')
        p('     GRS          0 h / 12 h / 24 h = %.2f (%.2f) / %.2f (%.2f) / %.2f (%.2f)\n', ...
          BD.GRSMean(1), BD.GRSSD(1), BD.GRSMean(2), BD.GRSSD(2), BD.GRSMean(3), BD.GRSSD(3));
        p('     GRS          global p = %.4f\n', BD.GRSP);
    end
    if isfield(BD, 'completionMean')
        p('     Steps done   0 h / 12 h / 24 h = %.2f (%.2f) / %.2f (%.2f) / %.2f (%.2f)  of 16\n', ...
          BD.completionMean(1), BD.completionSD(1), BD.completionMean(2), BD.completionSD(2), ...
          BD.completionMean(3), BD.completionSD(3));
        p('     Steps done   global p = %.4f; 0-12 p = %.4f, 0-24 p = %.4f, 12-24 p = %.4f\n', ...
          BD.completionP, BD.completionP12, BD.completionP13, BD.completionP23);
        p('     Last step s  %.0f (%.0f) / %.0f (%.0f) / %.0f (%.0f), global p = %.4f\n', ...
          BD.lastStepMean(1), BD.lastStepSD(1), BD.lastStepMean(2), BD.lastStepSD(2), ...
          BD.lastStepMean(3), BD.lastStepSD(3), BD.lastStepP);
    end
    if isfield(BD, 'timingRepairs') && ~isempty(BD.timingRepairs)
        p('\n   Timing sheet cells repaired or treated as not reached: %d\n', numel(BD.timingRepairs));
        for k = 1:numel(BD.timingRepairs)
            p('     %s\n', BD.timingRepairs{k});
        end
    end
catch ME_sec
    p('\n   [section could not be written: %s]\n', ME_sec.message);
    warning('rev01_write_report:section', ...
            'Section 4 participants and behaviour failed: %s', ME_sec.message);
end

p('\n%s\nEnd of report.\n', repmat('=', 1, 78));

end

%% ------------------------------------------------------------------------
function s = catLine(cats, cnts)
parts = arrayfun(@(i) sprintf('%s = %d (%.1f%%)', cats{i}, cnts(i), ...
                              100*cnts(i)/sum(cnts)), 1:numel(cats), ...
                 'UniformOutput', false);
s = strjoin(parts, ', ');
end

%% ------------------------------------------------------------------------
function printCell(fid, C, indent)
if isempty(C), fprintf(fid, '%s(unavailable)\n', indent); return; end
for i = 1:size(C, 1)
    fprintf(fid, '%s', indent);
    for j = 1:size(C, 2)
        v = C{i, j};
        if ischar(v) || isstring(v)
            fprintf(fid, '%-18s', char(v));
        elseif isnumeric(v) && isscalar(v)
            fprintf(fid, '%-18.5g', v);
        else
            fprintf(fid, '%-18s', '?');
        end
    end
    fprintf(fid, '\n');
end
end

% EOF
