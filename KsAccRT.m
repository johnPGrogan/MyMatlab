function h = KsAccRT(rt, acc, invertErrors, normaliseHeight, fillColours, ksArgs, plotArgs)
% function KsAccRT(rt, acc, invertErrors, fillColours, ksArgs, plotArgs)
% Plot KsDensity of the RT distributions, with correct/errors shown separately
% 
% Inputs:
%   rt  = [n x 1] vector of RTs
%   acc = [n x 1] vector of accuracies (1/0/NaN) - remove any ==2
%   invertErrors = 1=errors are mirrored in y-axis, 0=errors and correct both above x-axis
%   normaliseHeight = 1=normalise pdf of error/correct to both peak at 1, 
%     2=normalise by proportion of corr/err (e.g. 75% acc will mean corr is 3x higher)
%     0=don't normalise(default=0)
%   fillColours = 0=just line plot (uses colour order), 1=fill using default colours, 
%     or can give [2x3] rbg colours [correct; error] or cell {2 x 1} of
%     colour letters. fill uses 'hold on' to do error/correct separately
%   ksArgs = cell array that is passed into ksdensity as ksArgs{:}, can control
%     the x-axis limits via XI or Support (see help ksdensity) as well as bandwidth etc
%   plotArgs = cell array passed to plot/fill call. You can control the colour order
%     by using set(gca,'ColorOrder'...), if easier. Colour is controlled
%     separately for fill (fillColours argument), but transparency
%     ('FaceAlpha') should be controlled here for fill.
% 
% Outputs:
%   h = handle returned by plot() or two calls to fill()
% 

%% check inputs

if ~exist('invertErrors', 'var') || isempty(invertErrors)
    invertErrors = 1;
end
if ~exist('normaliseHeight', 'var') || isempty(normaliseHeight)
    normaliseHeight = 0;
end
if ~exist('fillColours', 'var') || isempty(fillColours)
    useFill = 0;
elseif isscalar(fillColours) % if was scalar, use that
    useFill = fillColours;
    if useFill % use default colours
        fillColours = get(gca,'ColorOrder');
    end
else
    useFill = 1; % use the colours from fillColours
end
if ~exist('ksArgs', 'var') || isempty(ksArgs)
    ksArgs = {};
end
if ~exist('plotArgs', 'var') || isempty(plotArgs)
    plotArgs = {};
end

%% plot it

for i = 1:2
    [y(i,:),x(i,:)] = ksdensity(rt(acc==(2-i)), ksArgs{:});
end

%% normalise?
if normaliseHeight==1 % make it peak at 1
    y = (y ./ max(y,[],2));

elseif normaliseHeight==2 % scale by p(correct)
    y = (y ./ max(y,[],2)) .* mean(acc == [1 0])';
end


%% invert?

if invertErrors
    y(2,:) = - y(2,:);
end

%% plot

if useFill
    x = [x, fliplr(x)]; % mirror it
    y = [y, zeros(size(y))]; % put zeros
    for i = 1:2
        h(i) = fill(x(i,:)', y(i,:)', fillColours(i,:), plotArgs{:});
        hold on;
    end

else % plot
    h = plot(x', y', plotArgs{:});

end







