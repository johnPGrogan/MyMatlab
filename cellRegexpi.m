function start = cellRegexpi(cellArray, expression)
% function x = cellregexpi(cellArray, expression)
% does regexpi on a cell or cell array
% converts the cell output into a matrix of indices found or zeros (if not
% found)
% 
% Input: 
%   cellArray = matrix of cells containing strings. no empty cells
%   pattern   = regexpi expression
% 
% Output:
%   start     = matrix of indices of start of expression in each cell, or 
%               zeros if not found. same size as cellArray
% 

if ~iscell(cellArray); error('cellArray is not a cell array'); end

start = regexpi(cellArray, expression); % run regexpi, gives indices or empty

start(cellfun(@isempty,start)) = {0}; % set empty to zero

start = cell2mat(start);