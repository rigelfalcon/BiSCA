function bshat = bshat_component(self, component, para_bs)

    if nargin < 3 || isempty(para_bs)
        para_bs = self.para_bs;
    end

    switch component
        case {'xi', 'xi_xi'}
            mask = false(size(self.para_bs.kernel.para.h));
            mask(1, :, :) = true;
            mask(:, 1, :) = true;
            para_bs.kernel.para.h(~mask) = 0;
        case {'alpha', 'alpha_alpha'}
            mask = true(size(self.para_bs.kernel.para.h));
            mask(1, 1, :) = false;
            para_bs.kernel.para.h(~mask) = 0;
        case {'xi+xi'}
            mask = false(size(self.para_bs.kernel.para.h));
            mask(1, 1, :) = true;
            para_bs.kernel.para.h(~mask) = 0;
        case {'alpha+alpha'}
            mask = true(size(self.para_bs.kernel.para.h));
            mask(1, :, :) = false;
            mask(:, 1, :) = false;
            para_bs.kernel.para.h(~mask) = 0;
        case {'xi+alpha'}
            mask = false(size(self.para_bs.kernel.para.h));
            mask(1, :, :) = true;
            mask(:, 1, :) = true;
            mask(1, 1, :) = false;
            para_bs.kernel.para.h(~mask) = 0;
        otherwise
            error('the component is not exsit, only xi and alpha');
    end

    bshat = pred_bs(self.f, para_bs, self.fxfy);

end