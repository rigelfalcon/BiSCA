function fit(self)
    self.bs = get_bs2fit(self);
    self.init_para();

    switch lower(self.optimization)
        case 'gmm'
            self.GMM();
        case 'nonlinearoptimization'
            self.NonlinearOptimization();
        otherwise
            error('no such method')
    end
    disp(['finish fit for name: [', num2str(self.name), ']']);

end