function I = nextOccurances(key,leads,varargin)
n = size(key,1);
reversedKey = key(n:-1:1);
reversedI = table.previousOccurances(reversedKey,leads,varargin{:});
I = reversedI(n:-1:1);
end