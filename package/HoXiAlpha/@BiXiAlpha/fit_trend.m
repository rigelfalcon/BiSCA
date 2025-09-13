function trend = fit_trend(self)
    % fit trend later use this trend to get the parameter of xi

    %% only fit xi ignore the peak, nonpara
    if strcmpi(self.method_xi, 'robust') % % robust detrend
        [~, ~, trend] = RobustDetrend(convert_to_double(mean(log(self.s), 2)), 3, 0.975, self.f); %LT= trenddecomp(self.s);

    elseif strcmpi(self.method_xi, 'linear') % % linear detrend
        n = 1;
        trend = detrend(log(mean(self.s, 2)), n);
    end

    % fit trend in log scale, addtive model in the raw scale
    trend = exp(trend(:));

end