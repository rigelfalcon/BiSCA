function varargout = convert_to_double(varargin)
% Convert input data to double if it is not already
varargout = cell(1, nargin);
for i = 1:nargin

    if isa(varargin{i}, 'double')
        varargout{i} = varargin{i};

    else
        if isa(varargin{i}, 'gpuArray')
            varargout{i} = gather(varargin{i}); % Transfer data from GPU to MATLAB workspace
        else
            varargout{i} = varargin{i};
        end
        varargout{i} = double(varargout{i});
    end
end
end