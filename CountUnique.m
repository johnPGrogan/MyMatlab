function [numInstances,uniqs] = CountUnique(mat,dim)
% function [numInstances,uniqs] = CountUnique(mat,dim)
% Count each instance of all instances

if nargin < 2
    dim = 'all';%default along rows
    dim2 = 1;
else 
    if length(dim)==1
        dim2=dim;
    else 
        dim2 = min(dim); % if multiple, take smallest
    end
end

uniqs = unique(mat);%get uniques
if iscell(mat)
    f = cellfun(@isnan, uniqs, 'UniformOutput',0);
    f = cellfun(@sum, f) > 0; 
    uniqs(f) = [];
else
    uniqs(isnan(uniqs)) = [];%remove NaNs;
end

% numInst = NaN(length(uniqs),size(mat,3-dim));%preallocate
numInstances = [];
if iscell(mat)%if cell
    
    for i = 1:length(uniqs)%count matches to each unique
        numInst = sum(strcmp(mat,uniqs(i)),dim);
        numInstances = cat(dim2,numInstances,numInst);
    end
    
elseif isnumeric(mat) || islogical(mat)%if number
    
    for i = 1:length(uniqs)
        numInst = sum(mat == uniqs(i),dim);
        numInstances = cat(dim2,numInstances,numInst);
    end
    
end