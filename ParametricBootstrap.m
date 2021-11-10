function [p, nullBetas, trueBeta] = ParametricBootstrap(dataTab, nullFormula, altFormula, nBoots, glmeArgs, tail, doPlot)
% function [p, nullBetas, trueBeta] = ParametricBootstrap(dataTab, nullFormula, altFormula, nBoots, glmeArgs, tail, doPlot)
% 
% Perform a parametric bootstrap (Bůžková, Lumley & Rice, 2011) in order to
% see how likely an interaction term is under a null model. We fit the null
% model (which is only missing the interaction term of interest) to the
% data, and then generate samples from this null model nPerms times,
% fitting the alternative model and keeping the coefficient estimate for
% the interaction term. This gives a distribution of the interaction from
% the null model which does not contain it.
% We then calculate the p-value as the probability of the true estimate
% given this null distribution.
% 
% Inputs:
%   dataTab = table to pass to fitglme for regression, each term in the
%       formulae must be a column. Variables should be z-scored
%   nullFormula = fitglme formula for the null model, which does not
%       contain the term of interest
%   altFormula = fitglme formula for the alternate model, which is
%       identical to the nullFormula except with one extra term
%   nBoots = number of bootstraps to do (default = 1000)
%   glmeArgs = cell array of any arguments to pass into fitglme (e.g.
%       {'link','logit','distribution','binomial'}. default is {}
%   tail = -1, 0 or 1. 0 = two-tailed (above or below zero), -1 is
%       one-tailed (below zero only) and 1 is one-tailed (above zero only).
%       default is 0.
%   doPlot = 0 or 1, 1 = plot hist of beta with true beta shown. default 0
% 
% Outputs:
%   p = p-value (probability of the alternative term-of-interest under the
%       null distribution
%   nullBetas = coefficient estimates of the alternate model
%     term-of-interest, generated from the null model [nBoots x 1]
%   trueBeta = alternate term-of-interest coefficient estimate from the
%     alterate model fit to the real data (i.e. the 'real' estimate)
% 
% John Grogan, 2021                                                                                     



%% check inputs
    
if ~istable(dataTab); error('dataTab must be a table'); end

if ~exist('nBoots','var') || isempty(nBoots); nBoots = 1000; end

if ~exist('glmeArgs','var') || isempty(glmeArgs)
    glmeArgs = {};
elseif ~iscell(glmeArgs)
    error('glmeArgs must be a cell or cell array');
end

if ~exist('tail','var') || isempty(tail)
    tail = 0; 
elseif ~any(tail == [-1 0 1])
    error('tail must be -1, 0 or 1');
end


%% do the base regressions and check they differ by one term only 

nullFit = fitlme(dataTab, nullFormula, glmeArgs{:});
altFit = fitlme(dataTab, altFormula, glmeArgs{:}); 

altInNull = ismember(altFit.CoefficientNames, nullFit.CoefficientNames);

if sum(~altInNull) ~= 1
    error('altFormula must contain ONE term that is not in the nullFormula');
end

altTermInd = find(~altInNull); % index of altTerm

% print term using
fprintf('running parametric bootstrap on the term %s, with %d bootstraps\n', altFit.CoefficientNames{altTermInd}, nBoots);

trueBeta = altFit.Coefficients.Estimate(altTermInd); % get true effect in data

%% do bootstrapping


regTab1 = dataTab; % copy this to overwrite data
nullBetas = zeros(nBoots,1); % preallocate

for i = 1:nBoots
    if mod(i,100)==0; disp(i); end % counter

    % generate samples from the null model
    regTab1.(nullFit.ResponseName) = random(nullFit);
    
    % fit altModel
    altFit1 = fitlme(regTab1, altFormula);

    % get beta
    nullBetas(i) = altFit1.Coefficients.Estimate(altTermInd);

end


%% do plot?

if doPlot
    figure();
    histogram(nullBetas);
    xline(trueBeta);
    xlabel('generated \beta from null model');
    ylabel('frequency');
end

%% calculate p-value


if tail == 1
    p = mean(nullBetas >= trueBeta); % prob of getting as big or bigger

elseif tail == -1
    
    p = mean(nullBetas <= trueBeta); % prob of getting as small or smaller

elseif tail == 0

    p = mean(abs(nullBetas) >= abs(trueBeta)); % prob of getting as extreme or more (i.e. pos and neg together)
    
end

end

