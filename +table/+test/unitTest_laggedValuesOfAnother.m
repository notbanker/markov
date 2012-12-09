function unitTest_laggedValuesOfAnother
ds = table.test.dsExample2(15);
lags = (1:2)';
ds.laggedTemperature = table.lagged.values(ds.temperature,lags,ds.dayOfWeekName);

f = @table.test.onMonTueWedReportTue;
[ds.strangeLaggedTemperature,I] = table.lagged.valuesOfAnother(ds.temperature,lags,f,ds.dayOfWeekName);

g = @table.test.onMonTueWedReportTue_Rain;
[ds.rainLaggedTemperature,I] = table.lagged.valuesOfAnother(ds.temperature,lags,g,ds.dayOfWeekName,ds.rain);

ds,
end