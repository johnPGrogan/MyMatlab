function myYlim(YData, lims)
% wrapper func to call ylim relative to minMax of h.YData or similar
% will adjust for negative values
% Inputs:
%   YData = vector of YData points (e.g. h.YData, [h.YData, h1.YData])
%   lims = relative limits e.g. [-0.05 0.05] to adjust the minMax by
% 

yl = minMax(YData,'all')'; %[min max]

ylim(yl .* (1 + lims .* sign(yl))); % flip sign of lims if yl is negative
