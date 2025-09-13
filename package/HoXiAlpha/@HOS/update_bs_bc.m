function update_bs_bc(self)
updateseg(self);

if self.update_auto && ~isempty(self.Xloc) && ~self.update_lock
    tic
    calc_bs_bc(self);
    self.update_auto = false;
    self.update_cross = true;
    t = toc;
    disp(['bs/bc calculation time: ', num2str(t), 's']);
end

end
