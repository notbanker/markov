function [rows,cols] = size(ds,k)
% function [rows,cols] = size(ds,k)

if isa(ds,'dataset'),
    [rows,cols] = size(ds);
elseif isa(ds,'struct'),
    vn = fieldnames(ds);
    rows = size(ds.(vn{1}),1);
    cols = length(vn);
else
    error('Expecting struct or dataset');
end

if nargin>1,
   switch k
       case 1 
           % rows = rows
       case 2
           rows = cols;
       otherwise
           error('Expecting dsSize(ds,1) or dsSize(ds,2)');
   end
end

end