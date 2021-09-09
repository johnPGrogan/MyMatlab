function x2 = tail(x, n)
if ~exist('n','var')
    n = 10;
end

dims = size(x);

nDims = length(dims);

nDimsBig = sum(dims > 1);

if nDimsBig==0
    x2 = x;
elseif nDimsBig < 3
    x2 = squeeze(x);
    n = min(n, size(x2,1));
    x2 = x2(end-n+1:end,:);
else
    [~,i] = sort(dims,'descend');
    x2 = permute(x,i);
    n = min(n, size(x2,1));
    x2 = x2(end-n+1:end,:);
    
end

