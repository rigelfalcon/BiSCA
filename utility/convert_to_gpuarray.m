function varargout = convert_to_gpuarray(varargin)
    % Converts multiple variables to gpuArray

    varargout = cell(1, nargin);
    force_cpu = strcmpi(getenv('BISCA_FORCE_CPU'), '1');
    use_gpu = ~force_cpu && (gpuDeviceCount("available") > 0);

    for i = 1:nargin
        varargout{i} = varargin{i};
    end

    if ~use_gpu
        return
    end

    try
        for i = 1:nargin
            % Check if the variable is already a gpuArray
            if isa(varargin{i}, 'gpuArray')
                varargout{i} = varargin{i}; % No need for conversion
            else
                varargout{i} = gpuArray(varargin{i}); % Convert to gpuArray
                if ispc
                    varargout{i} = single(varargout{i});
                end
            end
        end
    catch
        % Fallback to CPU arrays if GPU conversion fails (e.g., unsupported GPU toolchain)
        for i = 1:nargin
            varargout{i} = varargin{i};
        end
    end

end
