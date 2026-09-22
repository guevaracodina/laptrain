function snirfAux2csv(fName)
% SNIRFAUX2CSV  Export AUX channels from a SNIRF file to CSV
%
%   snirfAux2csv(fName)
%
%   INPUT
%       fName : string or char
%           Full path to .snirf file
%
%   OUTPUT
%       A .csv file saved in the same folder as fName, with:
%           - Column 1: time
%           - Columns 2–14: AUX channels (names from output.aux.name)
%
%   REQUIREMENTS
%       SnirfLoad must be in path

check_homer_path;
% --- Load SNIRF file ---
output = SnirfLoad(fName);

% --- Basic checks ---
try
    aux = output.aux;
catch
    error('No AUX field found in file.');
end

if isempty(aux)
    error('AUX field exists but is empty.');
end

nAux = numel(output.aux);
if nAux ~= 13
    warning('Expected 13 AUX channels, found %d.', nAux);
end

% --- Extract time vector (assumed common across AUX) ---
time = output.aux(1).time;

% --- Preallocate data matrix ---
nSamples = numel(time);
dataMat  = zeros(nSamples, nAux + 1);

% First column: time
dataMat(:,1) = time(:);

% --- Extract AUX data ---
varNames = cell(1, nAux + 1);
varNames{1} = 'time';

for idx = 1:nAux
    dataMat(:, idx+1) = output.aux(idx).dataTimeSeries(:);
    varNames{idx+1}   = output.aux(idx).name;
end

% --- Convert to table ---
T = array2table(dataMat, 'VariableNames', varNames);

% --- Define output CSV name ---
[filePath, fileBase, ~] = fileparts(fName);
% outName = fullfile(filePath, [fileBase '.csv']);

% --- Get parent folder name ---
[~, parentFolder] = fileparts(filePath);

% --- Construct new filename ---
outName = fullfile(filePath, [parentFolder '_' fileBase '.csv']);

% --- Write CSV ---
writetable(T, outName);

end