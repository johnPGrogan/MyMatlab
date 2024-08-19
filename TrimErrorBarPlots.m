function h = TrimErrorBarPlots(h, timesToTrim, reverseDirection)
% function h = TrimErrorBarPlots(h, times, timesToTrim, reverseDirection)
% 
% given a set of times (1 per line/patch), deleting XData,YData,Vertices,
% Faces from h where times(xdata) is > timesToTrim.
% can switch to remove before times.
% It will give a load of warnings during the editing, as the
% x/y/vertices/faces will not match up between edits.
%  
% Inputs:
%   h = errorBarPlot output {n, 2} with [lines patches] as columns, and n
%       different lines. Uses h.XData as times.
%   timesToTrim = [n, 1] vector of times to trim at
%   reverseDirection : 1=remove timepoints before timesToTrim, 0=remove
%     after (default = 0)
% 
% Outputs:
%  h = updated cell array of handles

if ~exist('reverseDirection','var') || isempty(reverseDirection)
    reverseDirection = 0;
end

% check sizes + dims
% timesToTrim  = [500; 700];

%% get times to keep

if reverseDirection
    trim = @(times, t) times > t; % return times above t
else
    trim = @(times, t) times < t; % return times below t
end


%% loop through lines/patches, deleting data outside them

n = size(h,1);

for i = 1:n

    for j = 1:2
%         tInds = trim(h{i,j}.XData,  timesToTrim(i));
        
        if j==2 % do vertices first
            h{i,j}.Vertices = h{i,j}.Vertices(trim(h{i,j}.XData,  timesToTrim(i)),:);
            h{i,j}.Faces = h{i,j}.Faces(trim(h{i,j}.XData,  timesToTrim(i)));
        end
        h{i,j}.YData = h{i,j}.YData(trim(h{i,j}.XData,  timesToTrim(i))); % then y
        h{i,j}.XData = h{i,j}.XData(trim(h{i,j}.XData,  timesToTrim(i))); % then x
    end
end





