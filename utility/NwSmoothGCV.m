function [yhat,nw_gcv]=NwSmoothGCV(x,yr,hList,xq,N)
[Tx,dx]=size(x);
[dh]=size(hList,1);
[Ty,dy]=size(yr);

if nargin<5|| isempty(N)
    N=ceil(1*(Tx^(1/dx)));
end
if nargin<4||isempty(xq)
    %     xlist=multispace(min(x),max(x),N,[],[],true);
    xq=get_ndgrid_scatter(x,'list',N);
end


% classUnderlying() may be unavailable without certain toolboxes; fall back to class().
try
    dtype = classUnderlying(x(1));
catch
    dtype = class(x(1));
end
[bsize] = bytesize(1,'b',dtype);

[~,freemem]=test_memory();
Nfree=freemem/bsize.b;
israndomize=size(x,1)^2*size(x,2)>Nfree/2;
Txr=floor(sqrt((Nfree/2)/size(x,2)));
if israndomize
    yhat=zeros(Txr,dy);
else
    yhat=zeros(Tx,dy);
end
pdof=zeros(dh,1);
gcvm=zeros(dh,1);
mse=zeros(dh,1);
df=zeros(dh,1);
% fpp=cell(dh);



for i=1:dh
    h=hList(i,:);
    if israndomize
        idx=randi(Tx,Txr,1);
        [yhat(:,i),L]=NwSmooth(x(idx,:),yr(idx,:),h,x(idx,:));
    else
        [yhat(:,i),L]=NwSmooth(x,yr,h,x);%,L
    end
    
    % fpp{i} = griddedInterpolant(multispace(min(x),max(x),N,[],'cell',true),reshape(yq,[N*ones(1,dx),1]),'spline','linear');
    % xcell=mat2cell(x,Tx,ones(1,dx));
    % yhat(:,i)=reshape(squeeze(fpp{i}(xcell{:})),[Tx,1]);
    pdof(i)=(1 - mean(diag(L)));
    df(i)=trace(L);
    %     gcvm(i)= Ty * sum((yr - yhat(:,i)) .^ 2) / ((Ty - df(i)) ^ 2);
    if israndomize
        mse(i,:)=mean((yr(idx,:) - yhat(:,i)).^2);
    else
        mse(i,:)=mean((yr - yhat(:,i)).^2);
    end
    gcvm(i)=mse(i,:)/(pdof(i).^ 2);
end
[~,idmin]=min(gcvm);
nw_gcv.idmin=idmin;
nw_gcv.gcvm=gcvm;
nw_gcv.pdof=pdof;
nw_gcv.df=df;
nw_gcv.mse=mse;
nw_gcv.hmin=hList(idmin,:);
nw_gcv.pdof_min=pdof(idmin);
% nw_gcv.fpp=fpp{i};
yhat=yhat(:,idmin);
% disp(['[nw_smooth]: [',num2str(dx),'d] bandwidth - hmin [',num2str(nw_gcv.hmin),'] is [',num2str(idmin),'th] of h'])
disp(['[nw_smooth]: [',num2str(dx),'d] regression' ])
disp([' - selected bandwidth - [',num2str(find(ismember(hList, nw_gcv.hmin.','rows'))),'th] h'])
disp([' - normalized scale h =',num2str(nw_gcv.hmin)])
disp([' - original scale h =',num2str(nw_gcv.hmin.*std(x))])
end








