function A = bmap(mA, idx, A)
    % compare object type of mA and A, if class(A(i)) is not the same as class(mA(i)) then cast it
    if ~strcmpi(class(mA), class(A))
        A = create_empty_objarray(size(A), mA);
    end

    A(idx) = mA;
end