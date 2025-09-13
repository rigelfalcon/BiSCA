
function [] = closematlabpool
if isempty(gcp('nocreate'))==0
    delete(gcp('nocreate'));
end