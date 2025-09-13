function init_xi(self)
    %% fit the trend
    %% init xi para
    init_para_xi(self);
    bound_xi(self);
    fit_xi(self);
end