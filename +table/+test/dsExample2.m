function ds = dsExample2(nSamples)
ds = dataset;
ds.temperature = ceil(50+rand(nSamples,1)*20);
ds.dayOfWeek = rem((1:nSamples)',7);
ds.rain = ceil(rand(nSamples,1)*2);
ds.dayOfWeekName = cell(nSamples,1);
for k=1:nSamples,
    ds.dayOfWeekName{k} = datestr(ds.dayOfWeek(k),'D');
end
end