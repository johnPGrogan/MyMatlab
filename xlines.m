function h = xlines(vals, varargin)
% plots multiple xlines - calls xline() on each value in val
% passes varargin to each

for i = 1:numel(vals)
    h(i) = xline(vals(i), varargin{:});
end
    