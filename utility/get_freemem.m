function bsize=get_freemem(onlyram)
    
if nargin<1 ||isempty(onlyram)
    onlyram=false;
end
    [~,freemem]=test_memory(0,onlyram);
    [bsize,scale] = bytesize( freemem, 'B','B');
    if nargout<1
        fprintf('Free memory: ');
        fprintf('%.2f %s\n', bsize.(scale), scale );
    end
end