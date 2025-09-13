% recover compacted data back to full matrix
function bs = get_bs2fit(self)

    switch self.compact
        case 'wedge'
            bs = wedge2full(self.bs, true);
        case 'tril'
            bs = tril2full(self.bs, 0, false);
        case 'quad1'
            bs = self.bs;
        case 'full'
            bs = self.bs;
    end

end