function unitTest_immutable

ds = dataset;
n = 50000;
ds.age = ceil(50*rand(n,1));
ds.height = 2+0.3*randn(n,1);
ds.pimples = ceil(5*rand(n,1));
ds.year_of_birth = ceil(1950+25*rand(n,1));
ds.year = 2012*ones(n,1);

% Supply all the fields needed
F = {'table.unitTests.determine_age','table.unitTests.hide_year_of_birth','table.unitTests.add_pimples'};
I = table.feval.equalPartition(5,ds.year_of_birth);
[dsOut,dsIn] = table.feval.immutable(F,ds,I);

% Supply less than what is needed
ds = table.rmfield(ds,'age');
[dsOut,dsIn] = table.feval.immutable(F,ds,I);




end