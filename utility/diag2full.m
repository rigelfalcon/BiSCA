function A=diag2full(D)
% take first dimenssion as the diag of frontal slice

sz=size(D);

A=zeros([sz(1),sz(1),prod(sz(2:end))]);
[p,~,n]=size(A);
idx=(1:p+1:p*p)'+((1:n)-1)*(p.^2);
A(idx)=D(:);
A=reshape(A,[sz(1),sz]);



end





