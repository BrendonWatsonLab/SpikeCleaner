function fileList = getAllExtFiles(dirName, fileExtension, includeSubfolders)
    % Function to get all files with a specific extension in a directory.
    % 
    % dirName: The directory to search in.
    % fileExtension: The file extension to look for (e.g., 'npy').
    % includeSubfolders: 1 to include subfolders, 0 to not include subfolders.

    if includeSubfolders
        dirData = dir(fullfile(dirName, ['**/*.' fileExtension]));  % Search including subfolders
    else
        dirData = dir(fullfile(dirName, ['*.' fileExtension]));  % Search in the specified folder only
    end

    fileList = {dirData.name}';

    % Prepend the directory to the file names
    for i = 1:length(fileList)
        fileList{i} = fullfile(dirData(i).folder, fileList{i});
    end
end
