function R = rightInverseSimple(L)
R = nan(nanmax(L),1);
for j=1:length(R),
    f = find(L==j,1,'first');
    if isempty(f),
        R(j) = nan;
    else
        R(j) = f(1);
    end
end

end