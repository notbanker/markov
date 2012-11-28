function ds = testFunc_add_pimples(ds)
  ds.pimples = rem(ds.age,5);
end