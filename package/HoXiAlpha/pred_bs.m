function [bshat] = pred_bs(f, para_bs, fxfy)

    para_bs.kernel.para.h(isnan(para_bs.kernel.para.h)) = 0;
    bshat = zeros(size(f, 1), size(f, 1));
    switch para_bs.type

        case 'tstudent'
            [bshat(para_bs.idx_fxfy)] = para_bs.kernel.eval(fxfy);
            bshat = tril2full(bshat, false, false);
        otherwise
            error('no such method')

    end

end