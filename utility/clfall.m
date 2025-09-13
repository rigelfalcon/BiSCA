function clfall
    % by Jan
    % https://ww2.mathworks.cn/matlabcentral/answers/478668-how-to-clear-not-close-all-the-opened-figures
    FigList = findall(groot, 'Type', 'figure');
    if ~isempty(FigList)
        for iFig = 1:numel(FigList)
            try
                clf(FigList(iFig));
            catch
            end
        end
    end
end