function [r, c] = GetSubPlotShape(n, dim)
% function [r, c] = GetSubPlotShape(n, dim)
% use ceil&sqrt to get good shape for subplots depending on number
% will be rectangular, roughly square
% 
% n = number of total plots
% dim = 1 or 2. 1 will make it longer than wide, 2 will be wider than logn
% (unless the numbers are equal obviously)

r = ceil(sqrt(n));
c = ceil(n / r);

if exist('dim','var')
    both = [r c];
    if dim==2 % make it wider
        r = min(both);
        c = max(both);
    elseif dim==1 % make it long
        r = max(both);
        c = min(both);
    else
        error('dim must be 1 or 2');
    end
end

if nargout == 1 % return both in one number
    r = [r c];
end
end