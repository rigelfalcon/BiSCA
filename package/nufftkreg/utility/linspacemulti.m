function x=linspacemulti(x_min,x_max,N)
m=max(max(numel(x_min),numel(x_max)),numel(N));
x_min=broadcast(x_min,m);
x_max=broadcast(x_max,m);
[N,flag_broad]=broadcast(N,m);
if flag_broad ||checksame(N)    
    x=ones(m,N(1),class(x_min));    
    for i=1:m
        x(i,:)=linspace(x_min(i),x_max(i),N(i));
    end
else
    x=cell(m,1);
    for i=m:-1:1
        x{i}=linspace(x_min(i),x_max(i),N(i));
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