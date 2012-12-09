function unitTest_leadingValuesOfAnother
ds = table.test.dsExample2(15);
lags = (1:2)';
ds.leadingTemperature = table.leading.values(ds.temperature,lags,ds.dayOfWeekName);

f = @table.test.onMonTueWedReportTue;
[ds.strangeLeadingTemperature,I] = table.leading.valuesOfAnother(ds.temperature,lags,f,ds.dayOfWeekName);

g = @table.test.onMonTueWedReportTue_Rain;
[ds.rainLeadingTemperature,I] = table.leading.valuesOfAnother(ds.temperature,lags,g,ds.dayOfWeekName,ds.rain);

ds,
end