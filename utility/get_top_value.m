function topValues = get_top_value(M, threshold)
    % GETTOPVALUES Calculates the top threshold value of each row in a matrix
    %   M: input matrix
    %   threshold: top threshold as a percentage of the number of columns
    %   topValues: output vector containing the top threshold value of each column

    % Get the number of rows in the matrix
    [numRows, numCols] = size(M);

    numNan = sum(isnan(M), 1);

    % Calculate the top threshold for each row
    if threshold ~= 1
        topThreshold = ceil((numRows - numNan) * (1 - threshold));
    else
        topThreshold = 1;
    end

    % Sort each column in descending order
    sortedM = sort(M, 1, 'descend');

    % Extract the top threshold value of each COLUMN
    idx = (topThreshold + numNan);
    idx = sub2ind([numRows, numCols], idx(:), (1:numCols)');
    topValues = sortedM(idx);

end