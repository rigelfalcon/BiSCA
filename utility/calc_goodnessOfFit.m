function rsquared = calc_goodnessOfFit(y, yfit)
    y = mean(y, 2, 'omitnan');

    % Calculate the goodness of fit
    SSE = sum((abs(y - yfit)).^2, "all", 'omitnan'); % Sum of Squared Errors
    SST = sum((abs(y - mean(y, 1, 'omitnan'))).^2, "all", 'omitnan'); % Total Sum of Squares
    rsquared = 1 - SSE / SST; % R-squared

end