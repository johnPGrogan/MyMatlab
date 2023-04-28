function [Y, I] = minMax(X, dim)
% [Y, I] = minMax(X, dim)
% return minimum and maximum values
% Inputs:   X = vector or matrix of values
%           DIM = dimension to operate along for min() and max(), can be
%               integer or 'all'.  if not given, will call min() and max() 
%               without dimension argument, which usually runs across the
%               first dimension (see 'help min')
% 
% Outputs:  Y = multidimensional array of min and max values, concatenated
%               along DIM dimension, i.e. dim=1: 2xN array; dim=2: Nx2
%               array
%           I = indices returned from min() and max(). will not work if  
%               DIM = 'all'
% 
% John Grogan

if ~exist('dim','var') || isempty(dim) % run min + max without dims argument
    yMin = min(X);
    yMax = max(X);
    if isscalar(yMin) && isscalar(yMax)
        Y = [yMin yMax];
    else
        Y = [yMin; yMax];
    end

else % if dim is given
    
    if isstr(dim) && strcmp(dim,'all') % cannot get inds for 'all'
        yMin = min(X, [], dim);    
        yMax = max(X, [], dim);
        I = [];
        
        % cat
        Y = cat(1, yMin, yMax);
    elseif length(dim) > 1 % need to use linear, means indices are across combined
        [yMin, iMin] = min(X, [], dim,'linear');
        [yMax, iMax] = max(X, [], dim,'linear');

        % cat across smallest dimension
        
        I = cat(min(dim), iMin, iMax);
        % cat
        Y = cat(min(dim), yMin, yMax);
    else
        [yMin, iMin] = min(X, [], dim);
        [yMax, iMax] = max(X, [], dim);
        I = cat(dim, iMin, iMax);
        
        % cat
        Y = cat(dim, yMin, yMax);
    end
    
end