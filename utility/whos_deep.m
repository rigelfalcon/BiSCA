
function sizeInMB = whos_deep(input, showDetails)
    if isempty(input)
        sizeInMB = 0;
        return;
    end

    warning('off', 'MATLAB:structOnObject');

    if isobject(input)&& ~ischar(input) && ~isstring(input)
        input = struct(input);
    end
    warning('on', 'MATLAB:structOnObject');

    if isstruct(input)
        fields = fieldnames(input);
        sizeInMB = 0;
        for i = 1:numel(fields)
            for j = 1:numel(input)
                sizeInMB = sizeInMB + whos_deep(input(j).(fields{i}), showDetails);
            end
        end
    elseif iscell(input)
        sizeInMB = 0;
        for i = 1:numel(input)
            sizeInMB = sizeInMB + whos_deep(input{i}, showDetails);
        end
    else
        info = whos('input');
        sizeInMB = info.bytes / 1024 / 1024;
    end
end
