
function x = symlog(x)
% x = sign(x) .* log10(abs(sign(x)*1+x)); %log
x = sign(x) .* log10(1+abs(x)); %log
end
