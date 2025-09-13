function X = tril2full(A,isvec,isconj)
% X = tril2full(mat2tril(ones(4),[],false))
if nargin<3||isempty(isconj)
    isconj=true;
end
if nargin<2||isempty(isvec)
    isvec=false;
end
if ~isvec
    [m,~,c] = size(A);
    idx=(1:m)+m*(0:m-1);
    idx=repmat(idx',[1,c]);
    idx=idx+(0:c-1)*m*m;
    X=zeros(m,m,c,class(A));
    for i=1:c
        if isconj
            X(:,:,i) = tril(A(:,:,i)) + tril(A(:,:,i))';
        else
            X(:,:,i) = tril(A(:,:,i)) + tril(A(:,:,i)).';
        end
    end
    X(idx) = A(idx);
else
    X=tril2full(vec2tril(A,0),false);
end

end
