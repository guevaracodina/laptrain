% Define source and target folders (modify paths as needed)
sourceFolder = fullfile('..','..','data','original');   % EDIT: raw Oxysoft export
targetFolder = fullfile('..','..','data','BIDSsource');   % EDIT if needed

for subNum = 1:30
    subName = sprintf('sub-%02d', subNum);
    for lapNum = 1:3
        lapFolderName = sprintf('lap0%d', lapNum);
        
        % Define the source SNIRF file for the current lap (e.g., 'Lap 1.snirf')
        sourceSNIRF = fullfile(sourceFolder, subName, sprintf('Lap %d.snirf', lapNum));
        
        % Destination directory: targetFolder\lap0X\sub-??\nirs
        destDir = fullfile(targetFolder, lapFolderName, subName, 'nirs');
        if ~exist(destDir, 'dir')
            mkdir(destDir);
        end
        
        % Copy and rename the SNIRF file: sub-??_task-lap0X_nirs.snirf
        destSNIRF = fullfile(destDir, sprintf('%s_task-lap0%d_nirs.snirf', subName, lapNum));
        copyfile(sourceSNIRF, destSNIRF);
        
        % Define possible source video file paths (either .mov or .mp4)
        sourceVideoMov = fullfile(sourceFolder, subName, sprintf('S%02d.%d.mov', subNum, lapNum));
        sourceVideoMp4 = fullfile(sourceFolder, subName, sprintf('S%02d.%d.mp4', subNum, lapNum));
        
        if exist(sourceVideoMov, 'file')
            videoSource = sourceVideoMov;
            videoExt = '.mov';
        elseif exist(sourceVideoMp4, 'file')
            videoSource = sourceVideoMp4;
            videoExt = '.mp4';
        else
            videoSource = '';
        end
        
        % Copy and rename the video file if it exists:
        % sub-??_task-lap0X_video.ext
        if ~isempty(videoSource)
            destVideo = fullfile(destDir, sprintf('%s_task-lap0%d_video%s', subName, lapNum, videoExt));
            if ~exist(destVideo, 'file')
                copyfile(videoSource, destVideo);
            else
                fprintf('%s already exists!\n',destVideo)
            end
        end
    end
end

%% Resting state
for subNum = 1:30
    subName = sprintf('sub-%02d', subNum);
    for restingNum = 1:3
        % Define the source SNIRF file for the current resting session
        sourceSNIRF = fullfile(sourceFolder, subName, sprintf('Resting %d.snirf', restingNum));
        
        % Define destination directory: targetFolder\resting0X\sub-??\nirs
        restingFolder = sprintf('resting0%d', restingNum);
        destDir = fullfile(targetFolder, restingFolder, subName, 'nirs');
        if ~exist(destDir, 'dir')
            mkdir(destDir);
        end
        
        % Define destination file name: sub-??_task-resting0X_nirs.snirf
        destSNIRF = fullfile(destDir, sprintf('%s_task-%s_nirs.snirf', subName, restingFolder));
        copyfile(sourceSNIRF, destSNIRF);
    end
end