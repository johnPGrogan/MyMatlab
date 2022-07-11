function [h, circles] = DrawCircles(coords, radius, useFill, varargin)
% [h, circles] = DrawCircles(coords, radius, useFill, varargin)
% plots a series of circles at the coordinates given, with radii given
% coords is either column of complex coordiantes i.e. complex(x, y), or an
% Nx2 matrix of [x,y] coordinates
% radius is a scalar radius or a vector of different radii for each circle
% useFill = 1, will fill circles. varargin will be passed to fill
% varargin gets passed into fill or plot
% 
% John Grogan, 2019

thetaRadians = (-180:1:180) * (pi/180);%convert degrees to radians

%get draw circles with radii
circleX = radius*sin(thetaRadians);
circleY = radius*cos(thetaRadians);

circle = complex(circleX, circleY); % make into complex

if ~isreal(coords)
    circles = circle + coords; % add to each complex coord of centre
    
elseif size(coords,2)==2% if a matrix of [x,y]
    
    coords = complex(coords(:,1), coords(:,2));
    circles = circle + coords; % add to each complex coord of centre
    
else
    
    error('coords should an Nx1 vector of complex coordinates, or an Nx2 matrix of [x,y] coords')
end

if exist('useFill','var') && useFill==1
    h = fill(real(circles)', imag(circles)', varargin{:});
else
    % permute and plot
    h = plot(permute(circles, [2,1]), varargin{:}); 
end
end