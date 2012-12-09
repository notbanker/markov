function Z = compose(varargin)
switch nargin,
    case 0,
        Z = [];
    case 1,
        Z = varargin{1};
    case 2,
        Z = composeTwo(varargin{1},varargin{2});
    otherwise
        tail = composeTwo(varargin{end-1},varargin{end});
        Z = table.indexing.compose(varargin{1:end-2},tail);
end

    function XY = composeTwo(X,Y)
        XY = nan(size(Y));
        legitY = ~isnan(Y) & Y>0;
        XY(legitY) = X(Y(legitY));
    end

end
