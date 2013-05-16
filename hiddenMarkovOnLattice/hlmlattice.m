function [lattice,m] = hlmlattice(y,varargin)
% Default lattice

if isempty(varargin),
    lattice = 200;
else
    lattice = varargin{1}; 
end
if length(lattice)==1,
    m = lattice;
    lb = min(y)-4*std(diff(y));
    ub = max(y)+4*std(diff(y));
    lattice = linspace(lb,ub,2*m+1);
end
m = length(lattice);
assert(rem(m,2)==1,'Lattice should be odd length'); 
m = (length(lattice)-1)/2;

end