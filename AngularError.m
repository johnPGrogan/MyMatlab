function error = AngularError(x,y,unit)
% function error = AngularError(x,y,unit)
% calculate circular angular error for two angles
% 
% Inputs:
%   x&y: matrix/vector/scalar angles in radians (-pi:pi) or degrees
%       (-180:180). should be same size or broadcastable
%   isRadians: 1 = is radians, 0 = isdegrees

if ~exist('unit','var') || isempty(unit)
    unit = 'radians';
end

switch unit
    case 'radians'
        % check limits
        if any(abs(x) > pi,'all') || any(abs(y) > pi,'all')
            warning('Your inputs seem to be in degrees, but you specified radians. please check');
        end
        
        % calculate
        error = mod(x - y + pi, 2*pi) - pi;
        
    case 'degrees'
        if all(abs(x) <= pi,'all') && all(abs(y) <= pi,'all')
            warning('Your inputs seem to be in radians, but you specified degrees. please check');
        end
        
        % calculate
        error = mod(x - y + 180, 360) - 180;
    otherwise
        error('units not recognised. Should be "radians" or "degrees"');
end


end