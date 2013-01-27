function ds = hide_year_of_birth(ds)
  table.rmfield(ds,'year_of_birth');
end