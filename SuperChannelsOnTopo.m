function h = SuperChannelsOnTopo(chanlocs, emarker2)
% function h = SuperChannelsOnTopo(chanlocs, emarker2)
% Superimpose channels onto a topoplot, using same inputs as topoplot()
% 
% chan_locs = an EEG.chanlocs structure (>> help readlocs or >> topoplot example)
% emarker2 = {markchans}|{markchans marker color size linewidth} cell array specifying 
%            an alternate marker for specified 'plotchans'. Ex: {[3 17],'s','g'} 
%            {default: none, or if {markchans} only are specified, then {markchans,'o','r',10,1}}
%   
%% defaults

EMARKER2CHANS = [];      % mark subset of electrode locations with small disks
EMARKER2 = 'o';          % mark subset of electrode locations with small disks
EMARKER2COLOR = 'r';     % mark subset of electrode locations with small disks
EMARKERSIZE2 = 10;      % default selected channel location marker size
EMARKER2LINEWIDTH = 1;

%% overwrite
EMARKER2CHANS = emarker2{1};
if length(emarker2) > 2
    EMARKER2 = emarker2{2};
end
if length(emarker2) > 3
    EMARKER2COLOR = emarker2{3};
end
if length(emarker2) > 4
    EMARKERSIZE2 = emarker2{4};
end
if length(emarker2) > 5
    EMARKER2LINEWIDTH = emarker2{5};
end

%% get pos to plot
[~, ~, Th, Rd] = readlocs( chanlocs );
Th = pi/180*Th;                              % convert degrees to radians

[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates

%% shrink to match 

rmax = 0.5;             % actual head radius - Don't change this!
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
headrad = rmax;  % (anatomically correct)

squeezefac = rmax/plotrad;
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
                          % to plot all inside the head cartoon
x    = x*squeezefac;    
y    = y*squeezefac;   


ELECTRODE_HEIGHT = 2.1;  % z value for plotting electrode information (above the surf)

hold on;
h = plot3(y(EMARKER2CHANS),x(EMARKER2CHANS),ones(size((EMARKER2CHANS)))*ELECTRODE_HEIGHT,...
       EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);