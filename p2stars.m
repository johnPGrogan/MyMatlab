function stars = p2stars(p, criteria, symbol)
% function stars = p2stars(p, criteria, symbol)
% converts p values into cell array of strings of stars. default criteria =
% * = p < .05
% ** = p < .01
% *** = p < .001
% **** = p < .0001
% 
% inputs: 
%   p = matrix of p values, any size
%   criteria = vector of decreasing p value criteria. default is [.05, .01, .001, .0001]
%   symbols = symbol to return, default is '*'
% 
% outputs:
%   stars = cell array matching size of p, with one '*' per criteria passed
% 


if ~exist('criteria','var') || isempty(criteria)
    criteria = [.05, .01, .001, .0001];
end
if ~exist('symbol','var') || isempty(symbol)
    symbol = '*';
else
    assert(ischar(symbol) && length(symbol)==1, 'symbol must be 1 character')
end

stars = cell(size(p)); % preset

for i = 1:length(criteria)
    stars(p <= criteria(i)) = {repmat(symbol, 1, i)};
end
