function X = toNumericArray(ds,varargin)
% Deal out tabular struct or dataset fields into array
fns = table.fieldnames(ds);
[~,obsLoc] = ismember(varargin(:),fns);
if all(obsLoc>0),
    n = length(varargin);
    X = zeros(size(ds,1),n);
    fat = false;
    for k=1:n,
        x_ = ds.(fns{obsLoc(k)});
        if isnumeric(x_),
            if size(x_,2)>1,
                fat = true;
                X = X(:,1:n-1);
            end
            if fat,
                X = [X,x_];
            else
                X(:,k) = x_;
            end
        else
            warning(['table.deal2array is ignoring non-numeric field ',fns{obsLoc(k)}]);
        end
    end
else
    warning(['table.deal2array was passed an invalid field']);
    f = find(obsLoc<1);
    varargin{f(:)},
    error('Given up');
end

end