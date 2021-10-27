function h = emptyLegend(n, opts, allOpts, names, legOpts)
% function h = emptyLegend(n, opts, cols, names)
% plots NaN so that you have plot handles for an empty plot which can be
% used for legends
% 
% Inputs:
%   n = number of lines/plots to make
%   opts = nx1 cell array, each containing one cell of arguments per line.
%   allOpts = cell array of options to be passed to all lines 
%   names = cell array of legend names
%   legOpts = cell array of options passed to legend
% 
% Example:
%   emptyLegend(2, {{'-r','LineWidth',2}, {'--','Color',[1 .2 .3]}}, {'Marker','x'}, {'test','control'}, {'Location','Best'});


    if ~exist('legOpts','var') || isempty(legOpts)
        legOpts = {};
    end
    if ~exist('allOpts','var') || isempty(allOpts)
        allOpts = {};
    end
    
    % plot each set of options as NaN
    for i = 1:n
        h(i) = plot(NaN, opts{i}{:}, allOpts{:});
        hold on;
    end
    
    legend(h, names(1:n), legOpts{:}); % make legend


end