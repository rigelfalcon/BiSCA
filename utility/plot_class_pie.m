function plot_class_pie(v, labelNames,opt)
    % plotClassPie: Creates a pie chart for the proportions of unique classes in vector v
    %
    % Inputs:
    %   v - A vector containing class labels (e.g., [1; 3; 4; 1; 3; 3; 4])
    %   labelNames (optional) - Cell array of class names corresponding to unique labels in v
    %                            (e.g., {'Class 1', 'Class 2', 'Class 3', 'Class 4'})
    %
    % Example:
    %   v = [1; 3; 4; 1; 3; 3; 4];
    %   plotClassPie(v, {'Class 1', 'Class 2 (Zero)', 'Class 3', 'Class 4'});

    % Check if labelNames is provided, otherwise use default names
    if nargin < 2
        labelNames = strcat('Class ', string(1:max(v)));
    end
    if nargin<3 || isempty(opt)
        opt=struct();
    end
    opt=set_defaults(opt,'cmp',[]);
    FontSize=8;

    % Define all possible class labels (assume labels are contiguous integers)
    uniqueClasses = 1:max(v);

    % Count occurrences of each class using histcounts
    counts = histcounts(v, [0.5:1:(max(uniqueClasses)+0.5)]);  % Counts for each class

    % Ensure the counts vector matches the labelNames length
    if length(labelNames) ~= length(counts)
        error('The number of labels does not match the number of classes in v');
    end
    
    % Prepare proportions for the pie chart (in percentages)
    total = sum(counts); % Total count
    proportions = counts / total * 100; % Convert to percentages
    
    % % Define colors for each class (soft colors for better readability)
    % colors = [0.4 0.6 1;     % Soft blue for Class 1
    %           0.8 0.8 0.8;   % Light gray for Class 2 (Zero)
    %           0.6 1 0.4;     % Soft green for Class 3
    %           1 0.6 0.4];    % Soft orange for Class 4
    if isempty(opt.cmp)
        cmp=flip(jet(length(labelNames)));
    else
        cmp=opt.cmp;
    end

    % Create pie chart
    h = pie(proportions, strcat(labelNames, ' (', string(counts), ')'));
    f=gcf;
    p = get(0, "MonitorPositions");
    % Position=p(2,:);
    % Position=[Position(1)+200,Position(2)+200,150,150];
    % Position=[200,200,150,150];
    % f.Position=Position;
    % Customize pie chart colors
    patchHandles = findobj(h, 'Type', 'Patch'); % Find pie slices
    for i = 1:length(patchHandles)
        patchHandles(i).FaceColor = cmp(i, :);
    end

    % Add title
    % title('Proportion of Classes', 'FontSize', FontSize,'Clipping','off');

    % Format the text labels on the pie chart for better readability
    textHandles = findobj(h, 'Type', 'Text');
    for i = 1:length(textHandles)
        % Add percentage values to the text labels
        percentText = sprintf('%.2f%%', proportions(i));
        textHandles(i).String = strcat(textHandles(i).String, '    ', percentText);  % Add percentage
        textHandles(i).FontSize = FontSize; % Adjust font size for readability
        textHandles(i).Clipping='off';
    end
end
