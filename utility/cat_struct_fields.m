function output = cat_struct_fields(input, field_names, istranspose, issingle)
% Initialize the output structure
output = struct();

if nargin < 2 || isempty(field_names)
    % Get the field names from the input structure
    field_names = fieldnames(input);
end

if nargin < 3 || isempty(istranspose)
    % Get the field names from the input structure
    istranspose = false;
end

if nargin < 4 || isempty(issingle)
    % Get the field names from the input structure
    issingle = false;
end

if isCharString(field_names)
    field_names = {field_names};
end

% Loop through each field
for i = 1:numel(field_names)
    % Get the name of the current field
    field_name = field_names{i};

    if contains(field_names, '.')
        withdot = true;
    else
        withdot = false;
    end

    if contains(field_names, '(')
        isfunction = true;
    else
        isfunction = false;
    end

    % Concatenate the data from the current field across all structures
    concatenated_data = [];

    for j = 1:numel(input)

        if ~withdot
            data = input(j).(field_name);
        else
            eval(['data = input(j).', field_name, ';']); %#ok<EVLDOT>
        end

        if istranspose
            data = data';
        end

        if isempty(data)
            continue;
        elseif isempty(concatenated_data)
            concatenated_data = data;
        else
            concatenated_data = cat(1, concatenated_data, data);
        end

    end

    % Store the concatenated data in the output structure
    if ~withdot
        output.(field_name) = concatenated_data;
    else
        eval(['output.', field_name, '=concatenated_data;']); %#ok<EVLDOT>
    end

end

if issingle && numel(field_names) == 1

    if ~withdot
        output = output.(field_name);
    else
        eval(['output=output.', field_name, ';']); %#ok<EVLDOT>
    end

end

end
