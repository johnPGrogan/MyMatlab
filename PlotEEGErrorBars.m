function h = PlotEEGErrorBars(data, varargin)
% function h = PlotEEGErrorBars(data, varargin)
% subplot each channel and plot EEG data (multiple bins) vs time
% Inputs:
%   data = eeg data [nPP nTimes nBins nChans]. will average over 1st dim
%       and plot one chan at a time
%   varargin = 'name', value pairs for optional parameters [default]
%       chanNames: cell array of names of channels, to pick ones to plot [numeric indices]
%       chansToPlot: cell array of chanNames to plot [all]
%       xTimes: [1xN] vector of sampling times [1:end]
%       timesToPlot: [start end] inclusive window of times to plot [1 end]
%       condNames: cell array of names of bins [empty]
%       yLabels: cell array of y-axis labels [{''}]
%       lineCols: [Nx3] matrix of line colours [RGB] to be set on each subplot [unchanged]
%       plotargs: cell of 'name',value pairs to pass into errorBarPlot [empty]
%       supertitle: title to put at top []
%       subplotInds: [row cols] for subplots. use 0 for no subplot. default is
%       [ceil(n/ceil(sqrt(n))),ceil(sqrt(n))].
%       xLabels: cell array of x-axis labels, either one per subplot, or
%           one repeated for all [{''}]
% 
% Outputs:
%   h = errorBarPlot figure handle


defaults = {'chanNames', cellstr(num2str([1:size(data,4)]')); % numbers
            'chansToPlot', []; % plot all
            'xTimes', 1:size(data,2); % 1:end
            'timesToPlot', [-Inf Inf]; % all 
            'condNames', {}; % nothing
            'yLabels', {''}; % nothing
            'lineCols', []; % default for figure
            'plotArgs', {}; % nothing
            'superTitle',[]; % nothing
            'subplotInds',[]; % none
            'yLine', 1; % yes
            'xLabels', {''}; % empty
            };
        
[chanNames, chansToPlot, xTimes, timesToPlot, condNames, yLabels, lineCols, plotArgs, superTitle, subplotInds, yLine, xLabels]...
    = parsepvpairs(defaults(:,1), defaults(:,2), varargin{:});
        
%%
if isempty(chansToPlot)
    chansToPlot = chanNames;
end
chanInds = find(ismember(chanNames, chansToPlot));
n = length(chanInds);
if isempty(subplotInds)
    subplotInds = GetSubPlotShape(n);
end

% if isempty(subplotInds) % make a rectangle
%     subplotInds = [ceil(n/ceil(sqrt(n))),ceil(sqrt(n))];
% end

timeInds = xTimes >= timesToPlot(1) & xTimes <= timesToPlot(2);
xTimes = xTimes(timeInds);
xlims = minMax(xTimes,2);
if numel(yLabels)==1 % repeat
    yLabels = repmat(yLabels, n,1);
end
if numel(xLabels)==1 % repeat
    xLabels = repmat(xLabels, n,1);
end
for i = 1:n
    if subplotInds ~= 0
        subplot(subplotInds(1), subplotInds(2),i)
    end
    if ~isempty(lineCols); set(gca,'ColorOrder',lineCols); end
    
    h = errorBarPlot(sq(data(:,timeInds,:,chanInds(i))),'area',1,'xaxisvalues',xTimes(1,:),'plotargs',plotArgs);
    hold on;
    xline(0,':k');
    if yLine
        yline(0,':k');
    end
    
    xlim(xlims);
    ylabel(yLabels{i}); 
    xlabel(xLabels{i}); 
    title(chanNames{chanInds(i)});

end

if ~isempty(condNames)
    legend([h{:,1}], condNames, 'Location','Best')
end

if ~isempty(superTitle)
    SuperTitle(superTitle);
end

end