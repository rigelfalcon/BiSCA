function [x] = nonlinear_fit(self, lostfun, x0, lb, ub, nonlcon)
    self.iiter = 0;
    tic

    x0 = strict_bound_x0(x0, lb, ub);
    [x0, lb, ub, self.para_fit.parainfo] = remove_fixed_para(x0, lb, ub);

    if isempty(x0)
        warning('no para to optimize');
        x = recover_fixed_para(x0, [], [], self.para_fit.parainfo);

    else

        A = [];
        b = [];
        Aeq = [];
        beq = [];
        HonorBounds = false;

        if nargin < 6 || isempty(nonlcon)
            nonlcon = [];
        end

        switch self.para_fit.optfun
            case 'fmincon'
                options = optimoptions('fmincon', 'MaxIterations', self.para_fit.MaxIterations, ...
                    'MaxFunctionEvaluations', self.para_fit.MaxFunctionEvaluations, 'Display', self.para_fit.Display, ...
                    'Algorithm', self.para_fit.Algorithm, 'UseParallel', self.para_fit.UseParallel, ...
                    'StepTolerance', self.para_fit.StepTolerance, ...
                    'FunctionTolerance', self.para_fit.FunctionTolerance, ...
                    'OptimalityTolerance', self.para_fit.OptimalityTolerance, 'HonorBounds', HonorBounds);
                self.model.loss_history = nan(self.para_fit.MaxIterations + 1000, 1);
                [x, loss] = fmincon(@(x) lostfun(x), ...
                    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            case 'lsqnonlin'
                options = optimoptions('lsqnonlin', 'MaxIterations', self.para_fit.MaxIterations, ...
                    'MaxFunctionEvaluations', self.para_fit.MaxFunctionEvaluations, 'Display', self.para_fit.Display, ...
                    'Algorithm', self.para_fit.Algorithm, ...
                    'UseParallel', self.para_fit.UseParallel, ...
                    'StepTolerance', self.para_fit.StepTolerance, ...
                    'FunctionTolerance', self.para_fit.FunctionTolerance, ...
                    'OptimalityTolerance', self.para_fit.OptimalityTolerance, ...
                    'FiniteDifferenceStepSize', self.para_fit.FiniteDifferenceStepSize);
                self.model.loss_history = nan(self.para_fit.MaxIterations + 1000, 1);

                if ~isempty(nonlcon)
                    [x, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(@(x) lostfun(x), ...
                        x0, lb, ub, A, b, Aeq, beq, nonlcon, options);
                else
                    [x, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(@(x) lostfun(x), ...
                        x0, lb, ub, [], [], [], [], [], options);
                end

            otherwise
                error('not finish other method');
        end

        %% compute the confidence interval
        Nx = sum(sum(abs(jacobian), 1) > eps, 2);

        [ci_compat, se] = nlpar_ci_se(x, residual, 'jacobian', jacobian);
        t = x ./ se;

        x = recover_fixed_para(x, [], [], self.para_fit.parainfo);

        parainfo = self.para_fit.parainfo;
        parainfo.x = NaN(size(parainfo.x));
        ci(:, 1) = recover_fixed_para(ci_compat(:, 1), [], [], parainfo);
        ci(:, 2) = recover_fixed_para(ci_compat(:, 2), [], [], parainfo);
        se = recover_fixed_para(se, [], [], parainfo);
        t = recover_fixed_para(t, [], [], parainfo);

        self.loss = resnorm;
        self.model.x = x;
        self.model.ci = ci;
        self.model.se = se;
        self.model.t = t;
        self.model.resnorm = resnorm;
        self.model.exitflag = exitflag;
        self.model.output = output;
        self.model.lambda = lambda;
        self.model.jacobian = sparse(jacobian);
        self.model.loss_history(isnan(self.model.loss_history)) = [];
        elapsedtime = toc;
        self.model.elapsedtime = elapsedtime;
        % display the parent calling mfilename
        name = parent_mfilename();
        idx_end = find(~isnan(self.model.loss_history), 1, 'last');
        disp(['nonlinear_fit [', name, '] Fitted ', num2str(Nx), ' para with: ', num2str(self.iiter), ' iters, elapsed time: ', num2str(elapsedtime), ' secs, loss: ', num2str(resnorm), ', step loss: ', num2str(self.model.loss_history(max(idx_end - 1, 1)) - self.model.loss_history(idx_end))]);
    end

end