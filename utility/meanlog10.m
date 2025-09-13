function mu=meanlog10(x,dim,varargin)
if nargin<2 || isempty(dim)
    dim=find(size(x)>1);
end

mu=10.^mean(log10(x),dim,varargin{:});


end

