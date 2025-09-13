function filenames = get_files_by_ext(folder, ext, maxDepth)
% get_files_by_ext: Return full filenames with a given extension, searching recursively.
%
% USAGE:
%   filenames = get_files_by_ext(folder, ext, maxDepth)
%
% INPUT:
%   - folder    : Path to the directory to start the search from.
%   - ext       : File extension to search for (e.g., '.m'). 
%                 If empty, all files will be returned.
%   - maxDepth  : Maximum folder depth to search recursively.
%                 - 0:   Only search the top-level 'folder'.
%                 - 1:   Search 'folder' and its immediate subdirectories.
%                 - Inf: Search with no depth limit.
%                 - Default is 0 if not provided.
%
% OUTPUT:
%   - filenames : Cell array of full paths to the found files.

% --- Input Parsing ---
% If maxDepth is not provided or is empty, set it to 0 (search current folder only)
if nargin < 3 || isempty(maxDepth)
    maxDepth = 0; 
end

% --- Find files in the current directory ---
if isempty(ext)
    % Get all items in the folder
    files = dir(folder);
    % Exclude special directories '.' and '..'
    files = files(~ismember({files.name}, {'.', '..'}));
    % Filter out all directories to keep only files
    files = files(~[files.isdir]);
else
    % Get all items matching the extension pattern
    files = dir(fullfile(folder, ['*' ext]));
    % It's possible for a directory to match the pattern (e.g., folder 'data.m'),
    % so we explicitly filter out directories from the results.
    files = files(~[files.isdir]);
end

% Generate full paths for the files found at the current level
filenames = fullfile({files.folder}', {files.name}');

% --- Recurse into subdirectories if depth allows ---
if (maxDepth > 0)
    % Get all subdirectories in the current folder
    subdirs = dir(folder);
    subdirs = subdirs([subdirs.isdir]);
    % Exclude special directories '.' and '..'
    subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

    % For each subdirectory, make a recursive call
    for i = 1:length(subdirs)
        currentSubDir = fullfile(folder, subdirs(i).name);
        % Recursively call the function with decremented depth
        filesInSubDir = get_files_by_ext(currentSubDir, ext, maxDepth - 1);
        % Append the newly found files to the main list
        filenames = [filenames; filesInSubDir];
    end
end

end
% function filenames = get_files_by_ext(folder, ext)
% %get_files_by_ext Return full filenames with given extension or all files.  
% %   FILENAMES = get_files_by_ext(FOLDER, EXT) returns the full filenames of 
% %   files in the given FOLDER that match the specified file extension EXT.
% %   If EXT is empty, return all files except . and ..  
% 
% if isempty(ext)
%     files = dir(folder);
%     files = files(~ismember({files.name},{'.','..'})); 
% else
%     files = dir(fullfile(folder, ['*' ext]));    
% end
% 
% filenames = fullfile({files.folder}', {files.name}');
% 
% end