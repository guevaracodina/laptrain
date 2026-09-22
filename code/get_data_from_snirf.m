function get_data_from_snirf(parentDir)
% GET_AUX_FROM_SNIRF  Recursively extract AUX data from all .snirf files
%
%   get_data_from_snirf(parentDir)
%
%   INPUT
%       parentDir : string or char
%           Root directory to search for .snirf files
%
%   DESCRIPTION
%       This function searches recursively within parentDir for all .snirf
%       files and applies snirf2csv to each file found. Output CSV files
%       are saved alongside the original .snirf files.
%
%   REQUIREMENTS
%       - snirf2csv.m must be in MATLAB path
%       - SnirfLoad must be available (e.g., Homer3)
%
%   EXAMPLE
%       get_data_from_snirf('C:\data\snirf_dataset')

    % --- Input validation ---
    if nargin < 1 || ~isfolder(parentDir)
        error('Input must be a valid directory.');
    end

    % --- Find all .snirf files recursively ---
    fileList = dir(fullfile(parentDir, '**', '*.snirf'));

    if isempty(fileList)
        warning('No .snirf files found in the specified directory.');
        return;
    end

    fprintf('Found %d .snirf files.\n', numel(fileList));

    % --- Process each file ---
    for k = 1:numel(fileList)
        try
            fName = fullfile(fileList(k).folder, fileList(k).name);
            fprintf('Processing (%d/%d): %s\n', k, numel(fileList), fName);

            snirf2csv(fName);

        catch ME
            fprintf('Failed: %s\nReason: %s\n', fileList(k).name, ME.message);
        end
    end

    fprintf('Processing complete.\n');

end