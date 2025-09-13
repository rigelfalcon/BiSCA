function folders = get_folders_in_directory(directory)
% Authors:
% - Ying Wang
% - Min Li
% - ChatGPT
% https://github.com/LMNonlinear
% https://github.com/rigelfalcon
% Date: April 11, 2022

% use the dir function to get a list of all files and folders in the directory
contents = dir(directory);

% initialize an empty list to store the folder names
folders = {};

% loop through each item in the contents list
for i = 1:length(contents)
    % check if the current item is a folder and exclude "." and ".." folders
    if contents(i).isdir && ~strcmp(contents(i).name,'.') && ~strcmp(contents(i).name,'..')
        % if it is a folder, add its name to the list
        folders{end+1} = contents(i).name; %#ok<AGROW>
    end
end

end
