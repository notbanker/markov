function unitTest_lagged_valuesOfAnother
ds = table.unitTests.dsExample2(15);
lags = (1:2)';
ds.laggedTemperature = table.lagged.values(ds.temperature,lags,ds.dayOfWeekName);

f = @table.unitTests.onMonTueWedReportTue;
[ds.strangeLaggedTemperature,I] = table.lagged.valuesOfAnother(ds.temperature,lags,f,ds.dayOfWeekName);

g = @table.unitTests.onMonTueWedReportTue_Rain;
[ds.rainLaggedTemperature,I] = table.lagged.valuesOfAnother(ds.temperature,lags,g,ds.dayOfWeekName,ds.rain);

ds,
end