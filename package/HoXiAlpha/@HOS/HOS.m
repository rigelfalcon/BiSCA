classdef HOS < handle & matlab.mixin.Copyable
    %%
    % The class inherits from the handle class, which allows the objects of this class to have properties that are passed by reference, i.e., modifications to these properties in one object will be reflected in other objects that reference them.
    % The class has several properties, including X (the input signal), Fs (the sampling frequency of the input signal), t (the time vector), N (the number of samples per segment), Ns (the number of segments), Ntrial (the number of trials), window (the window size), nfft (the FFT size), frange (the frequency range of interest), compact (a flag indicating whether to compute compact or full bicoherence), method (the method for spectral analysis), normalize (a flag indicating whether to normalize the spectra), normalization (the type of normalization to use), NW (the number of windows for computing averages), overlap (the overlap between consecutive windows), isseg (a flag indicating whether to treat the input signal as a single segment or multiple segments).
    % The class has several methods, including init_para (to initialize the parameters), set and get methods for accessing and modifying the properties, and a show_bc method for displaying the bicoherence of the input signal.
    % The updateval function is called by the set methods of some properties to update the value of the property and set the update flag to true, indicating that other properties that depend on the updated property need to be recomputed.
    % The update_bs_bc method is called by the get methods of some properties to update the bispectrum and bicoherence of the input signal, which are computed using the bispec function from the Signal Processing Toolbox.
    % The fx and fy properties are computed using the ndgrid function to create a meshgrid of frequencies for plotting the bicoherence.
    % Authors:
    % Ying Wang
    % Min Li
    % Pedro A. Valdes Sosa
    properties
        X
        Xloc
        Fs
        t
        Nt
        N
        Ns
        Ntrial
        Nvar
        window
        nfft
        frange
        compact
        method
        normalize
        normalization = 'haubrich';
        taper = 'SineScaled'; % 'Slepian' 'Sine' 'SineScaled'
        droplast = true
        NW
        overlap
        random_sample = false;
        num_random_segments
        tf_return_res = false;
        tf_return_seg = false;
        tf_return_var = true;
        log_var_s = true;
        info

        f
        w

        fx
        fy

        name

    end

    properties
        sloc
        svarloc
        sresloc
        Sloc
        Cloc
        coefloc
        coefrawloc
        bsloc
        bcloc
        bsvarloc
        bsresloc
        bsdenomloc
        W2loc
        W3loc
        floc
        wloc
        update_auto = true;
        update_cross = true;
        update_lock = false;
    end

    properties %(Dependent)
        s
        S
        C
        coef
        coefraw
        bs
        bc
        bsdenom
        W2
        W3
        svar
        bsvar
        sres
        bsres
    end

    methods

        function self = HOS(X, Fs)
            % X dim: [self.N, self.Ns, self.Ntrial, self.Nvar]
            if nargin < 1
                X = [];
            end

            if nargin < 2
                Fs = 1000;
            end

            self.X = X;
            self.Fs = Fs;
            init_para(self)
        end

        function init_para(self)
            self.t = [0:1 / self.Fs:(self.N - 1) / self.Fs];

            self.window = 256;
            self.nfft = 2 ^ nextpow2(self.window);
            self.frange = [0, 50];
            self.compact = true;
            self.method = 'pmtm';
            self.normalize = true;
            self.normalization = 'hagihira';
            self.NW = 2.5; %0.625
            self.overlap = 0.75;

        end

        %% setter

        function set.X(self, v)
            [self.X, self.update_auto] = updateval(self.X, v, self.update_auto); %#ok<*PROPLC>
            self.X = self.X - mean(self.X, 1);

        end

        function set.window(self, v)
            [self.window, self.update_auto] = updateval(self.window, v, self.update_auto);

        end

        function set.f(self, v)
            self.f = v;

        end

        function set.N(self, v)
            self.N = v;
        end

        function set.Ns(self, v)
            self.Ns = v;
        end

        function set.update_auto(self, v)
            self.update_auto = v;
        end

        function set.nfft(self, v)
            [self.nfft, self.update_auto] = updateval(self.nfft, v, self.update_auto); %#ok<*MCSUP>
        end

        function set.frange(self, v)
            [self.frange, self.update_auto] = updateval(self.frange, v, self.update_auto);
        end

        function set.compact(self, v)
            [self.compact, self.update_auto] = updateval(self.compact, v, self.update_auto);
        end

        function set.method(self, v)
            [self.method, self.update_auto] = updateval(self.method, v, self.update_auto);
        end

        function set.normalize(self, v)
            [self.normalize, self.update_auto] = updateval(self.normalize, v, self.update_auto);
        end

        function set.normalization(self, v)
            [self.normalization, self.update_auto] = updateval(self.normalization, v, self.update_auto);
        end

        function set.NW(self, v)
            [self.NW, self.update_auto] = updateval(self.NW, v, self.update_auto);
        end

        function set.overlap(self, v)
            [self.overlap, self.update_auto] = updateval(self.overlap, v, self.update_auto);
        end

        function set.tf_return_res(self, v)
            [self.tf_return_res, self.update_auto] = updateval(self.tf_return_res, v, self.update_auto);
        end

        function set.tf_return_seg(self, v)
            [self.tf_return_seg, self.update_auto] = updateval(self.tf_return_seg, v, self.update_auto);
        end

        function set.log_var_s(self, v)
            [self.log_var_s, self.update_auto] = updateval(self.log_var_s, v, self.update_auto);
        end

        %% getter

        function s = get.s(self)

            if ~self.update_lock
                update_bs_bc(self);
                s = self.sloc;
            else
                s = self.s;
            end

        end

        function svar = get.svar(self)

            if ~self.update_lock
                update_bs_bc(self);
                svar = self.svarloc;
            else
                svar = self.svar;
            end

        end

        function S = get.S(self)

            if ~self.update_lock
                update_bs_bc(self);
                update_S_C(self);
                S = self.Sloc;
            else
                S = self.S;
            end

        end

        function C = get.C(self)

            if ~self.update_lock
                update_bs_bc(self);
                update_S_C(self);
                C = self.Cloc;
            else
                C = self.C;
            end

        end

        function coef = get.coef(self)

            if ~self.update_lock
                update_bs_bc(self);
                coef = self.coefloc;
            else
                coef = self.coef;
            end

        end

        function coefraw = get.coefraw(self)

            if ~self.update_lock
                update_bs_bc(self);
                coefraw = self.coefrawloc;
            else
                coefraw = self.coefraw;
            end

        end

        function bc = get.bc(self)

            if ~self.update_lock
                update_bs_bc(self);
                bc = self.bcloc;
            else
                bc = self.bc;
            end

        end

        function bsvar = get.bsvar(self)

            if ~self.update_lock
                update_bs_bc(self);
                bsvar = self.bsvarloc;
            else
                bsvar = self.bsvar;
            end

        end

        function bsdenom = get.bsdenom(self)

            if ~self.update_lock
                update_bs_bc(self);
                bsdenom = self.bsdenomloc;
            else
                bsdenom = self.bsdenom;
            end

        end

        function bs = get.bs(self)

            if ~self.update_lock
                update_bs_bc(self);
                bs = self.bsloc;
            else
                bs = self.bs;
            end

        end

        function W2 = get.W2(self)

            if ~self.update_lock
                update_bs_bc(self);
                W2 = self.W2loc;
            else
                W2 = self.W2;
            end

        end

        function W3 = get.W3(self)

            if ~self.update_lock
                update_bs_bc(self);
                W3 = self.W3loc;
            else
                W3 = self.W3;
            end

        end

        function f = get.f(self)

            if ~self.update_lock
                update_bs_bc(self)
                f = self.floc;
            else
                f = self.f;
            end

            [self.fx, self.fy] = ndgrid(f);
        end

        function w = get.w(self)

            if ~self.update_lock
                update_bs_bc(self)
                w = self.wloc;
            else
                w = self.w;
            end

        end

        function sres = get.sres(self)

            if ~self.update_lock
                update_bs_bc(self)
                sres = self.sresloc;
            else
                sres = self.sres;
            end

        end

        function bsres = get.bsres(self)

            if ~self.update_lock
                update_bs_bc(self)
                bsres = self.bsresloc;
            else
                bsres = self.bsres;
            end

        end

        function set.update_lock(self, val)

            if val && ~self.update_lock
                fields = fieldnames(self);
                idx_loc = endsWith(fields, 'loc');
                fields_loc = fields(idx_loc);
                fields_used = erase(fields_loc, 'loc');

                for i = 1:length(fields_used)
                    self.(fields_used{i}) = self.(fields_loc{i});
                end

            end

            self.update_lock = val;
        end

        %%
        function runall(self)
            update_auto = self.update_auto; %#ok<*PROP>
            self.update_auto = true;
            self.f;
            update_bs_bc(self);
            self.update_auto = update_auto;
        end

        function compact_result(self)
            self.update_lock = true;
            self.X = [];
            self.Xloc = [];
            self.coef = [];
            self.coefloc = [];
            fields = fieldnames(self);
            idx_loc = endsWith(fields, 'loc');
            set_field(self, fields(idx_loc), []);

        end

        %%
        function show_s_mag(self, ivar)

            if nargin < 2 || isempty(ivar)
                ivar = 1;
            end

            if isnumeric(ivar)
                s = squeeze(self.s);
                plot(self.f(2:end), 10 .* log10(s(2:end, ivar)), '-o', 'LineWidth', 2) %./opt.nfft %db
            elseif strcmpi(ivar, 'mean')
                plot(self.f(2:end), 10 .* log10((self.s(2:end, :))), '-o', 'LineWidth', 2) %./opt.nfft %db
            elseif strcmpi(ivar, 'diag')
                plot(self.f(2:end), 10 .* log10(mean(self.s(2:end, :), 2)), '-o', 'LineWidth', 2) %./opt.nfft %db
            end

            title(['spectrum'])
        end

        function show_bc_mag(self, ivar)

            if nargin < 2 || isempty(ivar)
                ivar = 1;
            end

            fig = gcf;
            if isnumeric(ivar)
                surfbc(self.fx(2:end, 2:end), self.fy(2:end, 2:end), abs(self.bc(2:end, 2:end, ivar)));
            elseif strcmpi(ivar, 'mean')
                surfbc(self.fx(2:end, 2:end), self.fy(2:end, 2:end), abs(mean(self.bc(2:end, 2:end, :), 3, 'omitnan')))
            elseif strcmpi(ivar, 'diag')
                bcdiag = abs(full2diag(self.bc));
                bcdiag(isnan(bcdiag)) = 0;
                plot(self.f, bcdiag, 'rx');
                hold on;
                plot(self.f, mean(bcdiag, 2), 'go');
                hold on;
            end

            title('magnitude part of data bicoherence')

        end

        function show_bc_phase(self, ivar)
            fig = gcf;
            set(fig, 'Position', [499, 357, 450, 321]);
            surfbs(self.fx, self.fy, angle(self.bc(:, :, ivar)))
            title('magnitude part of data bicoherence')
            cmin = min(angle(self.bc), [], 'all');
            cmax = max(angle(self.bc), [], 'all');

            if cmin < cmax
                clim([cmin, cmax])
            end

        end

        function show_bs_mag(self, ivar)
            if isnumeric(ivar)
                surfbs(self.fx, self.fy, 10 .* log10(abs(self.bs(:, :, ivar))))
            elseif strcmpi(ivar, 'mean')
                surfbs(self.fx, self.fy, 10 .* log10(abs(mean(self.bs, 3))))
            elseif strcmpi(ivar, 'diag')
                bsdiag = 10 .* log10(abs(full2diag(self.bs)));
                bsdiag(isnan(bsdiag)) = 0;
                plot(self.f, bsdiag, 'rx');
                hold on;
                bsdenomdiag = 10 .* log10(abs(full2diag(self.bsdenom)));
                bsdenomdiag(isnan(bsdenomdiag)) = 0;
                plot(self.f, bsdenomdiag, 'bo')
            end

            title('magnitude bispectrum of data')
            sgtitle('spectrum and bispectrum of the data')
        end

        function show_bs_mag_times_phase(self, ivar)
            figure
            fig = gcf;
            set(fig, 'Position', [499, 357, 450, 321]);
            v = abs(self.bc) .* angle(self.bc(:, :, ivar));
            surfbs(self.fx, self.fy, v)
            title('magnitude part of data')
            sgtitle('bispectrum of the data')
            cmin = min(v, [], 'all');
            cmax = max(v, [], 'all');

            if cmin < cmax
                clim([cmin, cmax])
            end

        end

        function show_series(self, ivar)

            if isempty(self.X)
                return;
            end

            if nargin < 2 || isempty(ivar)
                ivar = 1;
            end

            if isnumeric(ivar)
                x = reshape(self.X(:, :, :, ivar), [], size(self.X(:, :, :, ivar), 4));
            elseif strcmpi(ivar, 'mean')
                x = reshape(mean(self.X(:, :, :), 4), [], size(self.X, 4));
            end

            Nt = size(x, 1);
            t = [0:1 / self.Fs:(Nt - 1) / self.Fs]';
            plot(t, x)
        end

        function show_embeding(self, lag, dim)

            if nargin < 2
                lag = [];
            end

            if nargin < 3
                dim = 3;
            end

            if isempty(lag)
                [~, lag] = phaseSpaceReconstruction(self.X, lag, dim);
            end

            timeIndices = 1:lag:length(self.X) - 2 * lag; % Time progression

            % Set up the colormap
            cmap = colormap(jet(length(timeIndices) - 1)); % Create a color map with enough colors

            hold on; % Keep all lines on the same figure

            if dim == 3
                % Plot segments for 3D trajectory with different colors
                for i = 1:(length(timeIndices) - 1)
                    x_segment = [self.X(timeIndices(i), 1), self.X(timeIndices(i + 1), 1)];
                    y_segment = [self.X(timeIndices(i) + lag, 1), self.X(timeIndices(i + 1) + lag, 1)];
                    z_segment = [self.X(timeIndices(i) + 2 * lag, 1), self.X(timeIndices(i + 1) + 2 * lag, 1)];
                    line(x_segment, y_segment, z_segment, 'Color', cmap(i, :), 'LineWidth', 2);
                end

                xlabel('x');
                ylabel('y');
                zlabel('z');
            else
                % Plot segments for 2D trajectory with different colors
                for i = 1:(length(timeIndices) - 1)
                    x_segment = [self.X(timeIndices(i), 1), self.X(timeIndices(i + 1), 1)];
                    y_segment = [self.X(timeIndices(i) + lag, 1), self.X(timeIndices(i + 1) + lag, 1)];
                    line(x_segment, y_segment, 'Color', cmap(i, :), 'LineWidth', 2);
                end

                xlabel('x');
                ylabel('y');
            end

            hold off;

            colorbar; % Show color scale
            colormap('jet'); % Use 'jet' colormap or choose another colormap

        end

        function show(self, ivar)
            if nargin < 2 || isempty(ivar)
                ivar = 'mean';
            end

            subplot(4, 1, 1)
            show_series(self, ivar)
            subplot(4, 1, 2)
            show_s_mag(self, ivar)
            subplot(4, 1, 3)
            show_bs_mag(self, ivar)
            subplot(4, 1, 4)
            show_bc_mag(self, ivar)
            set(gcf, 'Position', [200, 0, 350, 1200]);
        end

    end % method

end