function newStruct = keepFields(oldStruct, fieldsToKeep)
% function newStruct = keepFields(oldStruct, fieldsToKeep)
% wrapper function to call rmfield() on the fieldnames(oldStruct) that are
% not in fieldsToKeep. may be slower than other methods available for large
% structures, but should work with struct arrays etc. Will give errors if
% some fieldsToKeep were not found in oldStruct, or if all fields are to be
% removed.
% 
% Inputs:
%   oldStruct = structure or structure array
%   fieldsToKeep = string or cell of strings of fields to keep which match
%       fieldnames(oldStruct)
% Outputs:
%   newStruct = structure (same size as oldStruct) with only the
%       fieldsToKeep remaining
% 
% see help rmfield 
% John Grogan, 2023

if ~isa(oldStruct,'struct')
    error('oldStruct is not a structure');
end

if ~isstring(fieldsToKeep) && ~ischar(fieldsToKeep) && ~iscell(fieldsToKeep) && ~iscellstr(fieldsToKeep)
    error('fieldsToKeep must be a string or cell-string');
end

% get current fields
oldFields = fieldnames(oldStruct);

% match
toKeep = ismember(oldFields, fieldsToKeep);

if sum(toKeep) ~= length(fieldsToKeep)
    error('some fieldsToKeep were not found in fieldnames(oldStruct)');
elseif all(toKeep)
    error('all fields will be removed');
end


% remove those fields
newStruct = rmfield(oldStruct, oldFields(~toKeep));
