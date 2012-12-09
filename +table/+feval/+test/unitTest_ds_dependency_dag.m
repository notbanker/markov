function ok = unitTest_ds_dependency_dag

ds = dataset;
n = 50000;
ds.age = ceil(50*rand(n,1));
ds.height = 2+0.3*randn(n,1);
ds.pimples = ceil(5*rand(n,1));
ds.year_of_birth = ceil(1950+100*rand(n,1));
ds.year = 2012*ones(n,1);

F = {'testFunc_determine_age','testFunc_hide_year_of_birth','testFunc_add_pimples'};
[after,dsMin,In,Out,Removed] = table.feval.dependencyDag(F,ds);
ok = true;

end