function h = ylines(vals, varargin)
% plots multiple ylines - calls yline() on each value in val
% passes varargin to each

for i = 1:numel(vals)
    h(i) = yline(vals(i), varargin{:});
end
    