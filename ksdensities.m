function [f, x, u, ksinfo] = ksdensities(yData,varargin)
% function [f, x, u, ksinfo] = ksdensities(yData,varargin)
% loop through a 2D matrix and pass each column to ksdensity, will either
% plot it (with hold on), or will return the outputs of ksdensity per
% column.
% Can pass gca (figure/axis handle) as first argument like ksdensity, and
% can pass in any additional arguments to ksdensity which will be called
% for each column.
% Inputs:
%   yData: 2D matrix, each column will be passed to ksdensity separately
%   varargin: any other arguments to pass to ksdensity
% 
% Outputs:
%   (these are the same as ksdensity, with one row per column in yData)
%   f = vector of density values
%   x = set of points f is evaluated at. plot(x',f') to plot them
%   u = bandwidths used
%   ksinfo = structure array of ksdensity 
% 
% John Grogan, 2021

%% is there a fig/ax handle passed?

% copied from ksdensity.m
[axarg,varargin] = axescheck(yData,varargin{:});
if isempty(axarg)
    axarg = {};
else
    axarg = {axarg};
end
if isempty(varargin)
    yData = [];
else
    yData = varargin{1};
    varargin(1) = [];
end

%%

nCols = size(yData,2);
for i = 1:nCols
    if all(isnan(yData(:,i))); continue; end % skip if all NaN
    if nargout % if requesting args, will not plot
        [f(i,:),x(i,:),u(i),ksinfo(i)] = ksdensity(yData(:,i),varargin{:});
    else
        hold on;
        ksdensity(yData(:,i),varargin{:});
    end
end