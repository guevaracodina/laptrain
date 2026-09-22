function snirf2csv(fName)
% SNIRF2CSV  Export SNIRF time series to CSV with measurement labels
%
%   snirf2csv(fName)
%
%   INPUT
%       fName : string or char
%           Full path to .snirf file
%
%   OUTPUT
%       A .csv file saved alongside the input file with:
%           - Column 1: time
%           - Columns 2..end: dataTimeSeries channels
%             named as SxxDxxWxxx (zero-padded)
%
%   REQUIREMENTS
%       SnirfLoad in path

    % --- Load file ---
    output = SnirfLoad(fName);

    % --- Extract core data ---
    try
        time = output.data.time(:);
        data = output.data.dataTimeSeries;
        meas = output.data.measurementList;
        wavelengths = output.probe.wavelengths;
    catch
        error('Invalid SNIRF structure: missing required fields.');
    end

    % --- Consistency check ---
    [nTime, nMeas] = size(data);
    if numel(time) ~= nTime
        error('Time vector length does not match dataTimeSeries.');
    end

    % --- Preallocate ---
    dataMat = [time, data];
    varNames = cell(1, nMeas + 1);
    varNames{1} = 'time';

    % --- Generate column names ---
    for idx = 1:nMeas
        s = meas(idx).sourceIndex;
        d = meas(idx).detectorIndex;
        wIdx = meas(idx).wavelengthIndex;

        % Resolve wavelength value
        if wIdx > numel(wavelengths)
            error('Wavelength index exceeds available wavelengths.');
        end
        wVal = wavelengths(wIdx);

        % Zero-padded formatting
        varNames{idx+1} = sprintf('S%02dD%02dW%03d', s, d, round(wVal));
    end

    % --- Create table ---
    T = array2table(dataMat, 'VariableNames', varNames);

    % --- Output filename: [parentFolder]_fileBase.csv ---
    [filePath, fileBase, ~] = fileparts(fName);
    [~, parentFolder] = fileparts(filePath);
    outName = fullfile(filePath, [parentFolder '_' fileBase '.csv']);

    % --- Write CSV ---
    writetable(T, outName);

end

% EOF