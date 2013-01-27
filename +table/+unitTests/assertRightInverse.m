function assertRightInverse(L,R,dom)
rangeL = L(dom);
assert(~any(isnan(R(rangeL))));
for y = rangeL(:)',
    assert(L(R(y))==y);
end
end