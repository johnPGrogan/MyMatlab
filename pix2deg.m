function [deg, ppd] = pix2deg(pix, res, length, distance)
% function [deg, ppd] = pix2deg(pix, res, length, distance)
% convert from pixels into visual degrees
% pix = pixels to conver
% res = resolution in pixels (vert or horiz) 
% length = length of screen in cm (vert or horiz - MUST MATCH RES dimension)
% distance = viewing distance in cm
% 
% returns deg, and pixelsPerDeg
% 
ppd = pi * (res) / atan(length/ distance/2) / 360;

deg = pix ./ ppd;