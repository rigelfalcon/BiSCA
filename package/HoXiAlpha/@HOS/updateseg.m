function updateseg(self)
    update = self.update_auto && ~self.update_lock;

    if ~isempty(self.X) && update ...
            && ~isempty(self.window) && self.window > 0 ...
            && ~isempty(self.overlap) && self.overlap > 0
        [N, Ns, Ntrial, Nvar] = size(self.X);
        self.N = N;
        self.Ns = Ns;
        self.Ntrial = Ntrial;
        self.Nvar = Nvar;
        self.Nt = N * Ns;

        X = self.X;
        X = reshape(X, [N * Ns, Ntrial, Nvar]);
        [N, Ntrial, Nvar] = size(X);
        window = self.window; %number for the length of window
        idx = 1:round(window * (1 - self.overlap)):(N - window + 1); % self.N->size(X,1)
        self.Ns = length(idx);
        idx = (idx(:) + (0:(window - 1))).'; % idx of each window
        idx = idx(:) + (N) * (0:(Ntrial - 1)); % for trials
        idx = idx(:) + (N * Ntrial) * (0:(Nvar - 1)); % for Nvars
        idx = reshape(idx, window, [], Nvar); %

        %--- Added: Random segment sampling ---%
        if ~isempty(self.random_sample)
            idx_flat = idx;

            % Determine number of segments to sample
            if ~isempty(self.num_random_segments)
                num_segments = min(self.num_random_segments, size(idx_flat, 2));
            else
                num_segments = size(idx_flat, 2); % Default to all
            end

            % Randomly select segments
            selected_segments = sort(randperm(size(idx_flat, 2), num_segments));
            idx = idx_flat(:, selected_segments, :);

            % Update metadata (collapse variables into 1)

            self.Ns = num_segments;
            self.Ntrial = 1;
        end

        %--- End of added code ---%

        self.Xloc = X(idx);
        self.N = size(self.Xloc, 1);

        self.t = [0:1 / self.Fs:(self.N - 1) / self.Fs]';
    end

end