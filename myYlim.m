function myYlim(YData, lims)
% wrapper func to call ylim relative to minMax of h.YData or similar
% will adjust for negative values
% Inputs:
%   YData = vector of YData points (e.g. h.YData, [h.YData, h1.YData])
%   lims = relative limits to adjust the minMax by (default is [-.05 .05],
%       i.e. +/- 5%)
% 

if ~exist('lims','var') || isempty(lims)
    lims = [-.05 .05];
end

yl = minMax(YData,'all')'; %[min max]

ylim(yl .* (1 + lims .* sign(yl))); % flip sign of lims if yl is negative
