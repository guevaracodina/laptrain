function channelCounts = qtnirs_extract_good_ch_scanwise(figFile)
%% Extract bar data from QC_Group_good-channels_scanwise.fig
% figFile = 'QC_Group_good-channels_scanwise.fig';

% Open figure invisibly
hFig = openfig(figFile, 'invisible');

% Find axes
ax = findall(hFig, 'Type', 'axes');

% Find bar objects within axes
hBar = findall(ax, 'Type', 'bar');

% Check that at least one bar object was found
if isempty(hBar)
    close(hFig);
    error('No bar object was found in the figure.');
end

% If there are multiple bar objects, use the first one
% (for a simple bar plot this is usually sufficient)
y = hBar(1).YData;

% Convert to column vector of size nChannels x 1
channelCounts = y(:);

% Optional sanity check
% fprintf('Extracted vector size: [%d %d]\n', size(channelCounts,1), size(channelCounts,2));
% disp(channelCounts);

% Close figure
close(hFig);
end