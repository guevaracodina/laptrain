dirName = fullfile('..','..','data','original');   % EDIT: raw Oxysoft export
cd(dirName);
% List all folders starting with 'S' followed by any two characters
folders = dir('S??');
for i = 1:length(folders)
    if folders(i).isdir
        oldName = folders(i).name;
        % Ensure the folder name has at least three characters and that the 2nd and 3rd are digits
        if length(oldName) >= 3 && all(isstrprop(oldName(2:3), 'digit'))
            number = str2double(oldName(2:3));
            % Check that the number is within the desired range (01 to 30)
            if ~isnan(number) && number >= 1 && number <= 30
                newName = sprintf('sub-%02d', number);
                movefile(oldName, newName);
            end
        end
    end
end
