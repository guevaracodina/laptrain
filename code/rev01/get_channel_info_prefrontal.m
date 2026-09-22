function T = get_channel_info_prefrontal(csvFile)
% GET_CHANNEL_INFO_PREFRONTAL  Authoritative channel table for the 24-channel
% prefrontal Brite montage, in Homer3 measurement-list order.
%
%   T = get_channel_info_prefrontal()
%   T = get_channel_info_prefrontal(csvFile)
%
%   Neurophotonics NPH-260108-1, revision 01.
%   Single source of truth for channel index, source, detector, MNI coordinate
%   and anatomical label. Written because three files previously disagreed:
%   get_channels_from_template listed the montage detector-major, the published
%   optodes_MNI_coordinates.csv listed it source-major, and the hard-coded
%   label lists in conn_mat_3_group_comparison.m and
%   fc_network_measures_script_2025.m placed the short channels differently
%   again.
%
%   The ordering below was read directly from output.dc.measurementList by
%   rev01_diagnose and is detector-major. Short-separation channels are
%   indices 3 (S03D01) and 24 (S10D08); their MNI coordinates are NaN.
%
%   INPUT
%       csvFile  optional path to optodes_MNI_coordinates.csv. If given, the
%                file is read and checked against this table, and any
%                disagreement raises an error.
%
%   OUTPUT
%       T  24-row table: Channel, Source, Detector, X, Y, Z, Label, Region

ch  = (1:24)';
src = [1 4 3 1 2 4 5 4 5 2 5 6 7 5 6 6 7 8 9 6 8 8 9 10]';
det = [1 1 1 2 2 2 2 3 3 4 4 4 4 5 5 6 6 6 6 7 7 8 8  8]';

X = [ 63  44 NaN  48  20  29  29  20  19   8  17   1  -2  14  -3 -19 -21 -28 -41 -10 -19 -36 -48 NaN]';
Y = [ 24  26 NaN  41  56  35  45  26  36  50  52  51  55  45  44  46  50  38  44  36  28  31  37 NaN]';
Z = [ 17  28 NaN  14  14  27  27  38  38  24  26  22  13  31  28  19  10  20  10  35  36  25  15 NaN]';

Label = { ...
    'Frontal_Inf_Tri_R'; 'Frontal_Inf_Tri_R'; 'SSC'; 'Frontal_Inf_Tri_R'; ...
    'Frontal_Sup_Medial_R'; 'Frontal_Mid_R'; 'Frontal_Mid_R'; 'Frontal_Sup_R'; ...
    'Frontal_Sup_R'; 'Frontal_Sup_Medial_R'; 'Frontal_Sup_R'; 'Frontal_Mid_L'; ...
    'Frontal_Sup_Medial_L'; 'Frontal_Sup_R'; 'Frontal_Sup_L'; 'Frontal_Mid_L'; ...
    'Frontal_Sup_Medial_L'; 'Frontal_Mid_L'; 'Frontal_Inf_Tri_L'; ...
    'Frontal_Sup_Medial_L'; 'Frontal_Mid_L'; 'Frontal_Inf_Tri_L'; ...
    'Frontal_Inf_Tri_L'; 'SSC'};

Region = { ...
    'Right Inferior Frontal Gyrus, Triangular Part'
    'Right Inferior Frontal Gyrus, Triangular Part'
    'Short-separation channel'
    'Right Inferior Frontal Gyrus, Triangular Part'
    'Right Superior Medial Frontal Gyrus'
    'Right Middle Frontal Gyrus'
    'Right Middle Frontal Gyrus'
    'Right Superior Frontal Gyrus'
    'Right Superior Frontal Gyrus'
    'Right Superior Medial Frontal Gyrus'
    'Right Superior Frontal Gyrus'
    'Left Middle Frontal Gyrus'
    'Left Superior Medial Frontal Gyrus'
    'Right Superior Frontal Gyrus'
    'Left Superior Frontal Gyrus'
    'Left Middle Frontal Gyrus'
    'Left Superior Medial Frontal Gyrus'
    'Left Middle Frontal Gyrus'
    'Left Inferior Frontal Gyrus, Triangular Part'
    'Left Superior Medial Frontal Gyrus'
    'Left Middle Frontal Gyrus'
    'Left Inferior Frontal Gyrus, Triangular Part'
    'Left Inferior Frontal Gyrus, Triangular Part'
    'Short-separation channel'};

T = table(ch, src, det, X, Y, Z, Label, Region, 'VariableNames', ...
          {'Channel', 'Source', 'Detector', 'X', 'Y', 'Z', 'Label', 'Region'});

% Consistency with get_channels_from_template
[chIdxL, chIdxR] = get_channels_from_template('prefrontal');
sscExpected = setdiff(1:24, [chIdxL, chIdxR]);
sscHere     = find(strcmp(Label, 'SSC'))';
assert(isequal(sscExpected, sscHere), ...
       'Short-channel indices disagree: template [%s] vs table [%s].', ...
       num2str(sscExpected), num2str(sscHere));

if nargin >= 1 && ~isempty(csvFile) && isfile(csvFile)
    C = readtable(csvFile, 'VariableNamingRule', 'preserve', 'HeaderLines', 1);
    assert(height(C) == 24, 'CSV has %d rows, expected 24.', height(C));
    assert(isequal(C.Source(:), src) && isequal(C.Detector(:), det), ...
           'CSV source-detector order does not match the measurement list.');
    fprintf('Channel table verified against %s\n', csvFile);
end

end

% EOF
