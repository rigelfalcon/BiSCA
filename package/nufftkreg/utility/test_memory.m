function [isallow, freemem] = test_memory(bitSize, onlyram)
    if nargin < 1
        bitSize = 4e9;
    end

    if nargin < 2 || isempty(onlyram)
        onlyram = true;
    end

    if ismac
        % Code to run on Mac platform
        [r, w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        memsize = stats(1);
        freemem = (stats(3)) * 1024;
    elseif isunix
        % Code to run on Linux platform
        [r, w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        memsize = stats(1) * 1024;
        freemem = (stats(3)) * 1024;
    elseif ispc
        % Code to run on Windows platform
        [userview, systemview] = memory;
        if onlyram
            freemem = systemview.PhysicalMemory.Available;
        else
            freemem = userview.MemAvailableAllArrays;
        end
    else
        disp('Platform not supported')
    end

    if freemem > bitSize
        isallow = true;
    else
        isallow = false;
    end
end