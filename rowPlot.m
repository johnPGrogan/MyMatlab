function h = rowPlot(data,varargin)
% h = rowPlot(data,varargin)
% plot 1st row against 2nd row
% varargin are passed into plot

if size(data,1) ~= 2
    error('data does not have 2 rows')
end

h = plot(data(1,:),data(2,:),varargin{:});

end