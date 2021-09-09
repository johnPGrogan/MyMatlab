function [pix, ppd] = deg2pix(deg, res, length, distance)
% function [pix, ppd] = deg2pix(deg, res, length, distance)e)
% convert from visual degrees into pixels
% deg = visual degrees
% res = resolution in pixels (vert or horiz) 
% length = length of screen in cm (vert or horiz - MUST MATCH RES dimension)
% distance = viewing distance in cm
% 
% returns pix, and pixelsPerDeg
% 
ppd = pi * (res) / atan(length/ distance/2) / 360;

pix = deg .* ppd;