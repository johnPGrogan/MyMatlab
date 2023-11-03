function h = SwapXYErrorBars(h)
% given the output of errorBarPlot, swap the x ans y axes for everything

if length(h)>1
    % call on each one
    h = arrayfun(@SwapXYErrorBars, h);
    return;
end

if iscell(h)
    h = h{1};
end

% get fields to swap
fn = fieldnames(h);
xFields = fn(cellRegexpi(fn, '^X')>0);
yFields = fn(cellRegexpi(fn, '^Y')>0);

% remove any not in both
[xMatch, yMatch] = ismember(cellfun(@(x) x(2:end), xFields,'Uni',0), cellfun(@(x) x(2:end), yFields,'Uni',0));

xFields = xFields(xMatch);
yFields = yFields(yMatch(xMatch));


% make a copy
h1 = h; % I think this is a pointer? seems to change when h changes

% swap
for i = 1:length(xFields)
    orig.(xFields{i}) = h1.(xFields{i});
    h.(xFields{i}) = h1.(yFields{i});
    h.(yFields{i}) = orig.(xFields{i});
end
