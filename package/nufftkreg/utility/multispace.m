function x=multispace(x_min,x_max,N,space,type,istranspose)

if nargin<4 || isempty(space)
    space=@linspace;
end
if nargin<5 || isempty(type)
    type='array';
end
if nargin<6 || isempty(istranspose)
    istranspose=false;
end
if strcmp(func2str(space),'logspace')
    x_min=log10(x_min);
    x_max=log10(x_max);
end


m=max(max(numel(x_min),numel(x_max)),numel(N));
x_min=broadcast(x_min,m);
x_max=broadcast(x_max,m);
[N,flag_broad]=broadcast(N,m);
if (flag_broad ||checksame(N)) && strcmp(type,'array')
    x=ones(m,N(1),class(x_min));
    for i=1:m
        x(i,:)=space(x_min(i),x_max(i),N(i));
    end
    if istranspose
        x = x.';
    end
else
    x=cell(m,1);
    for i=1:m
        x{i}=space(x_min(i),x_max(i),N(i));
        if istranspose
            x{i}= x{i}.';
        end
    end
end

end
function [x,flag_broad]=broadcast(x,m)
flag_broad=isscalar(x);
if flag_broad
    x=repmat(x,[m,1]);
end
end

function  flag_same=checksame(x)
flag_same=~sum(abs(diff(x)));
end