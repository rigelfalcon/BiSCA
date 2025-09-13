function grid=get_ndgrid_scatter(x,type,N,space)
if nargin<2 || isempty(type)
    type='cell';
end
if nargin<3 || isempty(N)
    N=size(x,1);
end
if nargin<4 || isempty(space)
    space=@linspace;
end
if iscell(x)
    x = cell2mat2(x(:)');
end

x=multispace(min(x,[],1),max(x,[],1),N,space).';
grid=get_ndgrid(x,type);
end


