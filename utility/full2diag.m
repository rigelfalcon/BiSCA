
function [D,idx]=full2diag(A,isvec)
if nargin<2 || isempty(isvec)
    isvec=true;
end

sz=size(A);
if sz(1)~=sz(2)
    error('frontal slice should be square');
end
A=reshape(A,[sz(1),sz(2),prod(sz(3:end))]);
[p,~,n]=size(A);
if ndims(A)<=3
    idx=(1:p+1:p*p)'+((1:n)-1)*(p.^2);
else
    error('refer get_patch')
end
D=A(idx);

if isvec
    D=reshape(D,[sz(1),sz(3:end),1]);
else
    D=reshape(D,[sz(1),prod(sz(3:end)),1]);
    for i=1:size(D,2)
        A(:,:,i)=diag(D(:,i));
    end
    D=reshape(A,[sz,1]);
end

end