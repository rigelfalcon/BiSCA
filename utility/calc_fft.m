function [Y,w,f,Psd]=calc_fft(X,singleSide,average,Fs,nfft,isdemean)



if nargin<3 || isempty(average)
    average=false;
end
if nargin<2 || isempty(singleSide)
    singleSide=true;
end

%% Nnode,Nt,Ntrial
[Nnode,Nt,Ntrial]=size(X);



if nargin<6 || isempty(isdemean)
    isdemean=false;
end
if nargin<5 || isempty(nfft)
    nfft=Nt;
end
if nargin<4 || isempty(Fs)
    Fs=nfft;
end
if isdemean
    X=X-mean(X,2);
end

%% reshape
X=permute(X,[1,3,2]);%numNodes,numTrials,numFreqs
X=reshape(X,[Nnode*Ntrial,Nt]);


%% FFT

Y=fft(X,nfft,2);


%% Single side
if singleSide
    if rem(nfft,2)==0
        numFreqs=nfft/2+1;
    else
        numFreqs=(nfft+1)/2;
    end
    Y=Y(:,1:numFreqs);
else
    Y=fftshift(Y,2);
    numFreqs=nfft;
end
%% PSD
if nargout>3
    Psd = Y.*conj(Y)./Nt;
    Psd=reshape(Psd,[Nnode,Ntrial,numFreqs]);
    Psd=permute(Psd,[1,3,2]);%numNodes,numFreqs,numTrials
    Psd=Psd/4;
end


%% reshape back
Y=reshape(Y,[Nnode,Ntrial,numFreqs]);
Y=permute(Y,[1,3,2]);%numNodes,numFreqs,numTrials
%% if average
if average
    Y=mean(Y,3);
    if nargout>3
        Psd=mean(Psd,3);
    end
end
%% output W
if singleSide && nargout>1

    w = [0:numFreqs-1]'/(numFreqs-1)*pi;
else

    if rem(nfft,2)==0
        w = [-nfft/2:(nfft/2-1)]'/nfft*pi;
    else
        w = [-(nfft-1)/2:(nfft-1)/2]'/(nfft-1)*pi;
    end
end
%% output F

if nargout>2
    if singleSide

        f = [0:numFreqs-1]'/(numFreqs-1)*Fs;
    else

        if rem(nfft,2)==0
            f = [-nfft/2:(nfft/2-1)]'/nfft*Fs;
        else
            f = [-(nfft-1)/2:(nfft-1)/2]'/(nfft-1)*Fs;
        end
    end
else
    f=[];
end




end