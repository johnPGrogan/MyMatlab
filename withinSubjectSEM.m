function sem = withinSubjectSEM(data)
% function sem = withinSubjectSEM(data)
% Apply the Cousineau (2005) within-subject error bar correction. This
% removes between-subject variation before calculating SEM, and so returns
% SEM for within-subject variation only. 
% (omits nans from mean+std)
% Inputs:
%  data = n-dimensional matrix with participants as rows, other dimensions
%    can be time/condition etc. must have >1 columns/pages etc containing 
%    non-nan data (i.e. must be between + within variation), otherwise
%    it cannot split between/within, so will just return between SEM with a
%    warning
% 
% Outputs:
%  sem = n-dimensional matrix of SEM, with 1 row (other dimensions will
%    match data
% 
% John Grogan, 2024 (based on errorBarPlot from matlib)

sz = size(data);
nd = length(sz);


% get pp-level means, plus overall mean
% this loop will allow any number of dimensions
overallMean = mean(data,1, 'omitnan'); % take mean over participants, [1, sz(2:end)]
ppMean = data; % store as data, will take means in loop, skipping pp [nPP, 1, sz(3:end)]
n = data; % copy here too
for i = 2:nd % also get over all dimensions
    
    overallMean = mean(overallMean, i, 'omitnan'); % will end up as scalar
    ppMean = mean(ppMean, i, 'omitnan'); % will end up as [nPP, 1]

    % also track number of non-nans per dimension
    n = sum(~isnan(n),i);

end

if all(n < 2) % if fewer than one condition per subject
    warning('only one condition per subject, so cannot split between/within variation, returning overall SEM')
else
    % adjust the data
    data = data - ppMean + overallMean;
end

% now take sem across pps, 
sem = std(data,[],1, 'omitnan') ./ sqrt(sz(1));

