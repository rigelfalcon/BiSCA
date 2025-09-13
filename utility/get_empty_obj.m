function A=get_empty_obj(a)
if length(a)>1
    A(2)=a(1);
else
    A(2)=a;
end
A=A(1);
end