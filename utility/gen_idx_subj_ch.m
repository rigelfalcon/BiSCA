function idx = gen_idx_subj_ch(Nchs)
    % Validate input
    if ~isvector(Nchs) || any(Nchs < 0)
        error('Invalid input. Nchs must be a vector of positive integers or 0.');
    end
    
    % Get the total number of subjects
    Nsubj = numel(Nchs);
    
    % Initialize a cell array to store indices for each subject
    idx_cell = cell(Nsubj, 1);
    
    % Generate indices for each subject separately
    for i = 1:Nsubj
        if Nchs==0
            continue;
        end
        idx_cell{i} = [i * ones(Nchs(i), 1), (1:Nchs(i))'];
    end
    
    % Concatenate indices from all subjects
    idx = cat(1, idx_cell{:});
end
% function idx = generate_channel_indices(Nchs)
%     % Validate input
%     if ~isvector(Nchs) || any(Nchs < 1)
%         error('Invalid input. Nchs must be a vector of positive integers.');
%     end
    
%     % Get the total number of subjects
%     Nsubj = numel(Nchs);
    
%     % Generate indices using meshgrid and reshape
%     [subject_idx, channel_idx] = meshgrid(1:Nsubj, 1:max(Nchs));
%     idx = [subject_idx(:), channel_idx(:)];
    
%     % Remove rows with subject index greater than the number of channels for that subject
%     idx = idx(idx(:, 2) <= Nchs(idx(:, 1)), :);
% end
