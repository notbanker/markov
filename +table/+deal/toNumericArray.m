function X = toNumericArray(ds,fieldsToInclude)
% Deal out tabular struct or dataset fields into array
fns = table.fieldnames(ds);
[~,obsLoc] = ismember(fieldsToInclude,fns);
if all(obsLoc>0),
    nInclude = length(fieldsToInclude);
    X = zeros(table.size(ds,1),nInclude);
    fat = false;
    for k=1:nInclude,
        x_ = ds.(fns{obsLoc(k)});
        if isnumeric(x_),
            if size(x_,2)>1,
                fat = true;
                X = X(:,1:nInclude-1);
            end
            if fat,
                X = [X,x_];
            else
                X(:,k) = x_;
            end
        else
            warning(['table.deal.toNumericArray is ignoring non-numeric field ',fns{obsLoc(k)}]);
        end
    end
else
    warning(['table.deal2array was passed an invalid field']);
    f = find(obsLoc<1);
    varargin{f(:)},
    error('Given up');
end

end