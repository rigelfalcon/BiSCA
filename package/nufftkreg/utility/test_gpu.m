function [isallow,freemem]=test_gpu(bitSize,count,isglobal)
    %% 
    % isglobal might cause error, here it assume gpu memory not changed much
    % 
    
    if nargin<1 || isempty(bitSize)
        bitSize=4e9;
    end
    if nargin<2||isempty(count)
        count=1;
    end
    if nargin<3||isempty(isglobal)
        isglobal=false;
    end
    
    
    if isglobal
        global usegpu freegpumem
        if isempty(usegpu) || isempty(freegpumem) 
            [freemem, isallow] = check_gpu(count, bitSize);
            usegpu=isallow;
            freegpumem=freemem;
        else
            freemem=freegpumem;
            if freemem>bitSize
                isallow=true;
            else            
                isallow=false;
            end
        end
    else
        [freemem, isallow] = check_gpu(count, bitSize);
    end
    
    
    end
    
    function [freemem, isallow] = check_gpu(count, bitSize)
    
    if gpuDeviceCount>0    
        try
            gpuInfo=gpuDevice(count);
        catch
            fprintf('gpuDevice() has error, skipped.\n');
            gpuInfo=gpuDevice(1);
        end    
        freemem=gpuInfo.AvailableMemory;
        if freemem>bitSize
            isallow=true;
        else
            isallow=false;
        end
    else
        freemem=0;
        isallow=false;
    end
    end