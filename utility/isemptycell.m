function emptybool=isemptycell(A)
    emptybool=(cellfun(@isempty,A));%find
end