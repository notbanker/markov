function unitTest_createKeyChoice
ds = table.test.dsExample2(7);
f = @table.test.monday;
[key,U,R,C] = table.indexing.createKeyChoice(f,ds.dayOfWeekName,ds.rain);
end