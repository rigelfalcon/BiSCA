function bs = compact_bs(bs, compact)
    if ischar(compact) || isstring(compact)
        switch compact
            case 'wedge'
                bs = mat2wedge(bs, false);
            case 'tril'
                bs = mat2tril(bs, 0, false);
            case 'quad1'
                bs = bs;
            case 'full'
                bs = bs;
            otherwise
                error('no such compact')
        end
    else
        bs = bs(compact);
    end

end