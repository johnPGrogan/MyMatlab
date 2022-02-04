function statTable = StatsTableFromLMECells(lmeCells, colNames, statName, keepIntercept, rowNames)
% function statTable = StatsTableFromLMECells(lmeCells, colNames, statName, keepIntercept, rowNames)
% extract coefficient values from a cell vector of LME/GLME outputs
% store as a table
% Inputs:
%   lmeCells = cell vector of LME/GLME outputs. will ignore empty cells
%   colNames = cell array of strings. names to give columns of table created. 
%   statName = name of statistic to get from Coefficients. default is 'pValue'
%   keepIntercept = whether to keep the intercept or not (default = 0)
%   rowNames = cell/string of rownames to use, must match number of effects
%       extracted
% 
% Outputs:
%   stats = table with row per term and column per cell
% 

%% check inputs

if ~iscell(lmeCells) || ~isvector(lmeCells)
    error('lmeCells must be a cell array of LME/GLME outputs'); 
end

if exist('colNames','var') && ~isempty(colNames)
    if numel(colNames) ~= numel(lmeCells) || ~isvector(colNames)
        error('colNames size does not match lmeCells size');
    end
else
    colNames = [];
end

if ~exist('statName','var') || isempty(statName)
    statName = 'pValue';
end

if ~exist('keepIntercept','var') || isempty(keepIntercept)
    keepIntercept = 0;
end



% remove empties
emptyCells = cellfun(@isempty, lmeCells);
lmeCells(emptyCells) = [];
colNames(emptyCells) = [];

% check all are LMe/GLME type
if any( ~cellfun(@(x) isa(x, 'LinearMixedModel') || isa(x, 'GeneralizedLinearMixedModel'), lmeCells))
    error('not all filled cells are LME/GLME types');
end


%% get stat

stat = zeros(length(lmeCells{1}.CoefficientNames), length(lmeCells));

for i = 1:length(lmeCells)
    stat(:,i) = lmeCells{i}.Coefficients.(statName);
end

% convert to table
statTable = array2table(stat,'RowNames',lmeCells{1}.CoefficientNames);

if ~isempty(colNames)
    statTable.Properties.VariableNames = colNames;
end

if ~keepIntercept
    statTable(1,:) = [];
end    

% change row names?
% pick row names
if exist('rowNames','var') && ~isempty(rowNames)
    statTable.Properties.RowNames = rowNames;
end

