function U=atriu(A,k)
if nargin<2 || isempty(k)
    k=0;
end
A=flip(A,2);
U=triu(A,k);
U=flip(U,2);
end