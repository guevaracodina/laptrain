for i = 4:30
    parentFolder = sprintf('sub-%02d', i);
    if ~exist(parentFolder, 'dir')
        mkdir(parentFolder);
    end
    mkdir(fullfile(parentFolder, 'nirs'));
end
