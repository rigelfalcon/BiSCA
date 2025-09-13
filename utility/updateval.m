function [v,update]=updateval(objv,v,update)

tf=isequal(objv,v);
if ~tf
    update=true;
else
    update=update|false;
end



end