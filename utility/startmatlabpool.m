
function [pool] = startmatlabpool(size,resources)
% Authors:
% Ying Wang

pool=[];
opened = 0;
poolinfo=gcp('nocreate');
if isempty(poolinfo)~=1
    opened = 1;
end
if nargin<2 ||isempty(resources)
    resources='local';
end
if nargin<1 || isempty(size)
    size=feature('numcores');
end


if ~opened %not opened
    if nargin==0
        pool=parpool(resources);
    else
        try
            pool=parpool(resources,size);%matlabpool('open','local',size);
        catch ce
            pool=parpool(resources);%matlabpool('open','local');
            size = pool.NumWorkers;
            disp(ce.message);
            disp(strcat('restart. wrong  size=',num2str(size)));
        end
    end    
else
    disp('matlabpool has started');
    if nargin==1 && size~=poolinfo.NumWorkers
        closematlabpool;
        startmatlabpool(size);
    else
        disp('continue old matlabpool');
    end    
end

