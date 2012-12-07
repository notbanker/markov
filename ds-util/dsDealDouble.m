function X = dsDealDouble(ds,varargin)
% Deal out struct or dataset fields into array
fns = dsNames(ds);
[~,obsLoc] = ismember(varargin(:),fns);
if all(obsLoc>0),
    n = length(varargin);
    X = zeros(size(ds,1),n);
    fat = false;
    for k=1:n,
        x_ = ds.(fns{obsLoc(k)});
        if size(x_,2)>1,
            fat = true;
            X = X(:,1:n-1);
        end
        if fat,
            X = [X,x_];
        else
            X(:,k) = x_;
        end
    end
else
    warning(['Invalid field name requested ']);
    f = find(obsLoc<1);
    varargin{f(:)},
    error('Given up');
end

end