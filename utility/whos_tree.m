function whos_tree(input, name, indent, showDetails, threshold)
if nargin < 5 || isempty(threshold)
    threshold = 1; % default threshold: 1 MB
end
if nargin < 4|| isempty(showDetails)
    showDetails = true;
end
if nargin < 3|| isempty(indent)
    indent = '';
end
if nargin < 2|| isempty(name)
    name = 'root';
end

sizeInMB = whos_deep(input, showDetails);
if sizeInMB >= threshold % only print members >= threshold MB
    fprintf('%s%s: %.2f MB\n', indent, name, sizeInMB);
end
warning('off', 'MATLAB:structOnObject');
if isobject(input)&& ~ischar(input) && ~isstring(input)
    input = struct(input);
end
warning('on', 'MATLAB:structOnObject');
if (isstruct(input) && isscalar(input)) || (isstruct(input)&&numel(input) >1 && showDetails)
    fields = fieldnames(input);
    for i = 1:numel(fields)
        for j = 1:numel(input)
            whos_tree(input(j).(fields{i}), fields{i}, [indent '  '], showDetails, threshold);
        end
    end
elseif (iscell(input) && isscalar(input)) || (iscell(input) && numel(input) >1 && showDetails)
    for i = 1:numel(input)
        whos_tree(input{i}, sprintf('cell{%d}', i), [indent '  '], showDetails, threshold);
    end
end
end
