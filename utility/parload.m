function dataCell = parload(inputFiles, uniformStruct, maxWorkers)
    % PARLOAD - Loads multiple files in parallel using parfor with memory management.
    %
    %   dataCell = parload(inputFiles, uniformStruct, maxWorkers)
    %
    % Inputs:
    %   inputFiles    - Cell array of file paths to load.
    %   uniformStruct - Optional boolean flag (default: false). If true,
    %                   attempts to return a structure array if all loaded
    %                   files contain structs with the same fields.
    %   maxWorkers    - Optional integer (default: automatic). Limits the
    %                   number of parallel workers to control memory usage.
    %
    % Outputs:
    %   dataCell - Cell array containing loaded data (structs or empty for
    %              missing files), or a structure array if uniformStruct is
    %              true and all structs have the same fields.
    %
    % Notes:
    %   - Loading many large files in parallel may consume significant memory.
    %     Use maxWorkers to limit simultaneous loading, or set to 1 for
    %     sequential loading if memory is a concern.
    %   - Parallel processing may use more memory than a sequential for loop
    %     due to data transfer overhead.
    %
    % Examples:
    %   files = {'file1.mat', 'file2.mat', 'file3.mat'};
    %   data = parload(files);            % Default parallel loading
    %   data = parload(files, true, 2);   % Limit to 2 workers, uniform output

    % Validate input
    if ~iscell(inputFiles)
        error('inputFiles must be a cell array.');
    end

    % Set defaults
    if nargin < 2
        uniformStruct = false;
    end
    if nargin < 3
        maxWorkers = []; % Use default pool size
    elseif ~isscalar(maxWorkers) || maxWorkers < 1
        error('maxWorkers must be a positive integer.');
    end

    % Start or adjust parallel pool
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(maxWorkers)
            parpool('Processes');
        else
            parpool('Processes', maxWorkers);
        end
    elseif ~isempty(maxWorkers) && pool.NumWorkers > maxWorkers
        delete(pool);
        parpool('Processes', maxWorkers);
    end

    % Initialize output
    dataCell = cell(1, numel(inputFiles));

    % Load files in parallel
    parfor i = 1:numel(inputFiles)
        if exist(inputFiles{i}, 'file')
            dataCell{i} = load(inputFiles{i});
        else
            warning('File %s not found. Skipping.', inputFiles{i});
            dataCell{i} = [];
        end
    end

    % Convert to structure array if requested
    if uniformStruct
        nonEmpty = ~cellfun(@isempty, dataCell);
        if any(nonEmpty)
            firstStructIdx = find(nonEmpty, 1);
            firstFields = fieldnames(dataCell{firstStructIdx});
            isUniform = true;
            for i = 1:numel(dataCell)
                if nonEmpty(i)
                    if ~isstruct(dataCell{i}) || ~isequal(fieldnames(dataCell{i}), firstFields)
                        isUniform = false;
                        break;
                    end
                end
            end
            if isUniform
                template = cell2struct(cell(size(firstFields)), firstFields);
                structArray = repmat(template, 1, numel(dataCell));
                for i = 1:numel(dataCell)
                    if nonEmpty(i)
                        structArray(i) = dataCell{i};
                    else
                        structArray(i) = template;
                    end
                end
                dataCell = structArray;
            else
                warning('Non-uniform structs detected. Returning cell array.');
            end
        else
            warning('All files missing or empty. Returning empty cell array.');
        end
    end

    % Optional: Uncomment to shut down pool and free resources
    % delete(gcp('nocreate'));
end
%{
function dataCell = parload(inputFiles, uniformStruct)
    % PARLOAD - Loads multiple files in parallel using parfor.
    %
    %   dataCell = parload(inputFiles, uniformStruct)
    %
    % Inputs:
    %   inputFiles    - Cell array of file paths to load.
    %   uniformStruct - Optional boolean flag (default: false). If true,
    %                   attempts to return a structure array if all loaded
    %                   files contain structs with the same fields.
    %
    % Outputs:
    %   dataCell - Cell array containing loaded data (structs or empty for
    %              missing files), or a structure array if uniformStruct is
    %              true and all structs have the same fields.
    %
    % Example:
    %   files = {'file1.mat', 'file2.mat', 'file3.mat'};
    %   data = parload(files); % Returns cell array
    %   data = parload(files, true); % Returns struct array if uniform

    % Validate input
    if ~iscell(inputFiles)
        error('inputFiles must be a cell array.');
    end

    % Set default for uniformStruct
    if nargin < 2
        uniformStruct = false;
    end

    % Start a parallel pool if not already active
    if isempty(gcp('nocreate'))
        parpool('Processes');
    end

    % Initialize output
    dataCell = cell(1, numel(inputFiles));

    % Load files in parallel
    parfor i = 1:numel(inputFiles)
        if exist(inputFiles{i}, 'file')
            dataCell{i} = load(inputFiles{i});
        else
            warning('File %s not found. Skipping.', inputFiles{i});
            dataCell{i} = [];
        end
    end

    % If uniformStruct is true, attempt to convert to structure array
    if uniformStruct
        % Check if all non-empty elements are structs
        nonEmpty = ~cellfun(@isempty, dataCell);
        if any(nonEmpty)
            % Get fields of first non-empty struct
            firstStructIdx = find(nonEmpty, 1);
            firstFields = fieldnames(dataCell{firstStructIdx});
            
            % Check if all non-empty elements are structs with same fields
            isUniform = true;
            for i = 1:numel(dataCell)
                if nonEmpty(i)
                    if ~isstruct(dataCell{i}) || ~isequal(fieldnames(dataCell{i}), firstFields)
                        isUniform = false;
                        break;
                    end
                end
            end

            if isUniform
                % Initialize structure array with empty template
                template = cell2struct(cell(size(firstFields)), firstFields);
                structArray = repmat(template, 1, numel(dataCell));
                
                % Populate structure array
                for i = 1:numel(dataCell)
                    if nonEmpty(i)
                        structArray(i) = dataCell{i};
                    else
                        structArray(i) = template; % Empty struct for missing files
                    end
                end
                dataCell = structArray; % Replace cell array with struct array
            else
                warning(['Not all files contain structs with the same fields. ' ...
                         'Returning cell array instead.']);
            end
        else
            warning('All files are missing or empty. Returning empty cell array.');
        end
    end

    % Optional: Shut down the pool (uncomment if desired)
    % delete(gcp('nocreate'));
end
%}
% function dataCell = parload(inputFiles)
%     % PARALLELLOADFILES - Loads multiple files in parallel using parfor.
%     %
%     %   dataCell = parallelLoadFiles(inputFiles)
%     %
%     % Inputs:
%     %   inputFiles - Cell array of file paths to load.
%     %
%     % Outputs:
%     %   dataCell - Cell array containing loaded data (structs or empty for missing files).
%     %
%     % Example:
%     %   files = {'file1.mat', 'file2.mat', 'file3.mat'};
%     %   data = parallelLoadFiles(files);
% 
%     % Validate input
%     if ~iscell(inputFiles)
%         error('inputFiles must be a cell array.');
%     end
% 
%     % Start a parallel pool if not already active
%     if isempty(gcp('nocreate'))
%         parpool('Processes');
%     end
% 
%     % Initialize output
%     dataCell = cell(1, numel(inputFiles));
% 
%     % Load files in parallel
%     parfor i = 1:numel(inputFiles)
%         if exist(inputFiles{i}, 'file')
%             dataCell{i,:} = load(inputFiles{i});
%         else
%             warning('File %s not found. Skipping.', inputFiles{i});
%             dataCell{i,:} = [];
%         end
%     end
% 
%     % Optional: Shut down the pool (uncomment if desired)
%     % delete(gcp('nocreate'));
% end