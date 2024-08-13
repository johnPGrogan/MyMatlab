function h = SuperChannelsOnTopoLine(chanlocs, chanInds, plotArgs, scalingFactor)
% function h = SuperChannelsOnTopoLine(chanlocs, emarker2)
% Superimpose channels onto a topoplot, using same inputs as topoplot()
% 
% chan_locs = an EEG.chanlocs structure (>> help readlocs or >> topoplot example)
% chanInds = indices of channels, will plot ellipse
% plotArgs = cell array of linestyles to pass into plot3 as plotArgs{:}
% scalingFactor = [optional] scaling of radius so you can expand/shrink (default=1, no scaling)
% 
% adapted from topoplot.m

if ~exist('plotArgs','var') || isempty(plotArgs)
    plotArgs = {};
end
if ~exist('scalingFactor','var') || isempty(scalingFactor)
    scalingFactor = 1;
end

%% get pos to plot
[~, ~, Th, Rd] = readlocs( chanlocs );
Th = pi/180*Th;                              % convert degrees to radians

[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates

%% shrink to match 

rmax = 0.5;             % actual head radius - Don't change this!
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
% headrad = rmax;  % (anatomically correct)

squeezefac = rmax/plotrad;
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
                          % to plot all inside the head cartoon
x    = x*squeezefac;    
y    = y*squeezefac;   

ELECTRODE_HEIGHT = 2.1;  % z value for plotting electrode information (above the surf)

%% make ellipse

a = diff(minMax(y(chanInds)))/2; % horizontal axis
b = diff(minMax(x(chanInds)))/2; % vertical axis

%% find centre

coords = complex(y(chanInds), x(chanInds));
[~, centreInds] = min(coords - mean(coords));
centre = coords(centreInds);

x0=real(centre); % x0,y0 ellipse centre coordinates
y0=imag(centre);
t=-pi:0.01:pi;
x1=x0+a*cos(t)*scalingFactor;
y1=y0+b*sin(t)*scalingFactor;

%% plot

hold on; % make sure on
h = plot3(x1, y1, repmat(ELECTRODE_HEIGHT,size(x1)),...
    plotArgs{:});
