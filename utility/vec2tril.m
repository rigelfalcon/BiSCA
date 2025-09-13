function Atril=vec2tril(A,ind,N)
% example: 3*10, from vec to get lower trangle 2*2*10
% Yingwang 5/10/2020
%  ind= 0 is the main diagonal
%  ind > 0 is above the main diagonal
%  ind < 0 is below the main diagonal.
if nargin<2 ||isempty(ind)
    ind=0;
end
[n,m]=size(A);  %
if nargin<3||isempty(N)
    if ind==0
        N=(-1+sqrt(1+8*n))/2;
    else
        error('need input size of square matrix N when the ind is not 0')
    end
end
% if ind == 0
%     l=(-1+sqrt(1+8*n))/2; % original formula
% else
    % l=(-1+sqrt(1+8*(n-abs(ind))))/2; % adjust for negative ind
% end
% l=(-1+sqrt(1+8*(n-abs(ind))))/2; % adjust for negative ind
onesind=tril(true(N),ind);
AA=zeros(N,class(A));
Atril=zeros(N,N,m,class(A));
for i=1:m
    AA(onesind)=A(:,i);
    Atril(:,:,i)=AA;
end

end

% AA=AA'
% n=(-1+sqrt(1+8*length(A)))/2 ;
% switch ind
% case -1
%     onesind=find(tril(ones(a),-1));
% case 1
%     onesind=find(tril(ones(a),1));
% case 0
%     onesind=find(tril(ones(a)));
% end