function group = gen_loopgroup(num_all, id_in_each, isuniform)
    if nargin < 3 || isempty(isuniform)
        isuniform = true;
    end
    num_groups = ceil(num_all / id_in_each);
    if isuniform
        id_in_each_floor = floor(num_all / num_groups);
        remainder = mod(num_all, num_groups);
        sizes = repmat(id_in_each_floor, num_groups, 1);
        if num_groups > 1
            sizes(1:remainder) = sizes(1:remainder) + 1;
        end
    else
        remainder = mod(num_all, id_in_each);
        sizes = repmat(id_in_each, num_groups, 1);
        if remainder ~= 0
            sizes(end) = remainder;
        end
    end
    group = mat2cell((1:num_all), 1, sizes)';
end