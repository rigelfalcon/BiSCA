function [bc,bs,f,w,psd]=call_compute_bispectrum(X,Fs,nfft,frange,compact,window,overlap)

segleng=window;
segshift=round(overlap*segleng);
epleng=size(X,1);
if (rem(nfft,2) == 0)
    maxfreqbins = nfft/2-1;
    wraw = [0:maxfreqbins-1]'/nfft*pi;
else
    maxfreqbins = (nfft-1)/2;
    wraw = [0:maxfreqbins-1]'/(nfft-1)*pi;
end


X=X(:);% each column is the data for a channel in this function 
[bs,bc,Pyy] = compute_bispectrum(X,segleng,segshift,epleng,maxfreqbins);
bs=permute(bs,[2,3,1]);
bc=bs./permute(bc,[2,3,1]);



fraw=wraw./pi.*Fs;
if ~isempty(frange)
    [fx,fy]=ndgrid(fraw);
    bcf=griddedInterpolant(fx,fy,bc,'nearest','none');
    bsf=griddedInterpolant(fx,fy,bs,'nearest','none');
    if ~compact && frange(1)>=0
        frange=[-frange(2),frange(2)];
    end
    fidx = find(fraw >= frange(1) & fraw <= frange(2));
    f=fraw(fidx);
    w=wraw(fidx);
    psd=Pyy(fidx);
    NFFT=length(f);
    [fx,fy]=ndgrid(f);
    bc=bcf(fx,fy);
    bs=bsf(fx,fy);
end



if compact
    [~,idx]=mat2tril(bc,0,false);
    bc(setdiff(1:NFFT*NFFT*1,find(idx)))=NaN;
    bs(setdiff(1:NFFT*NFFT*1,find(idx)))=NaN;
end

end










