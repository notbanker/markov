function R = rightInverse(L,dom)
rangeL = L(dom);
R = nan(nanmax(rangeL),1);
isInDom = ismember((1:length(L))',dom);
for j=1:length(R),
    f = find(L==j & isInDom,1,'first');
    if isempty(f),
        R(j) = nan;
    else
        R(j) = f(1);
    end
end

end