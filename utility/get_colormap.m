function cm=get_colormap(v)
if sum(v>0,"all")&&sum(v<0,"all")
    cm=turbo;
else
    cm=(hot);%flip
end

end