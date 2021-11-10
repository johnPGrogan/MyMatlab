function [p, trueBeta, origP, nullBetas] = ParametricBootstrap(dataTab, nullFormula, altFormula, nBoots, glmArgs, tail, doPlot, fitFunc)
% function [p, trueBeta, origP, nullBetas] = ParametricBootstrap(dataTab, nullFormula, altFormula, nBoots, glmArgs, tail, doPlot fitFunc)
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
%
% Inputs:
%   dataTab = table to pass to fitglme for regression, each term in the
%       formulae must be a column. Variables should be z-scored
%   nullFormula = fitglme formula for the null model, which does not
%       contain the term of interest
%   altFormula = this can be either a formula that is the same as
%       nullFormula except with one extra fixed effect term, OR it can be
%       the extra fixed effect term which will be added to the nullFormula.
%       e.g. 'y ~ 1 + x1 + x2 + x1:x2' or 'x1:x2'
%       identical to the nullFormula except with one extra fixed term
%   nBoots = number of bootstraps to do (default = 1000)
%   glmArgs = cell array of any arguments to pass into fitting function, e.g.
%       {'link','logit','distribution','binomial'}. default is {}
%   tail = -1, 0 or 1. 0 = two-tailed (above or below zero), -1 is
%       one-tailed (below zero only) and 1 is one-tailed (above zero only).
%       default is 0.
%   doPlot = 0 or 1, 1 = plot hist of beta with true beta shown. default 0
%   fitFun = which function to use for fitting the data: @fitlm, @fitlme,
%       @fitglm, or @fitglme. They take increasing times (in the order
%       written there) so best to use the simplest one you can, depending
%       on the model formula and glmArgs given. Default is to pick the
%       quickest one that will work (i.e. fitlm if no random effects or
%       glmArgs, but fitglme if random effects and glmArgs). If you are
%       passing in glmArgs but want to use fitglm, then specify it.
%
% Outputs:
%   p = p-value (probability of the alternative term-of-interest under the
%       null distribution
%   trueBeta = alternate term-of-interest coefficient estimate from the
%     alterate model fit to the real data (i.e. the 'real' estimate)
%   origP = alternate term-of-interest p-value from fitting the alternate
%       model to the read data (i.e. 'real' p-value)
%   nullBetas = coefficient estimates of the alternate model
%     term-of-interest, generated from the null model [nBoots x 1]
%
% John Grogan, 2021



    %% check inputs
    
    if ~istable(dataTab); error('dataTab must be a table'); end
    
    if ~exist('nBoots','var') || isempty(nBoots); nBoots = 1000; end
    
    if ~exist('glmArgs','var') || isempty(glmArgs)
        glmArgs = {};
    elseif ~iscell(glmArgs)
        error('glmArgs must be a cell or cell array');
    end
    
    if ~exist('tail','var') || isempty(tail)
        tail = 0;
    elseif ~any(tail == [-1 0 1])
        error('tail must be -1, 0 or 1');
    end
    
    if isempty(regexp(altFormula, '~','ONCE')) % if it is not a full formula, just one term
        % then make it into a full formula
        % need to insert it between ' ~ 1' and and random effects
    
        parenInd = regexp(nullFormula, '(','ONCE'); % get first parenthesis
    
        if isempty(parenInd) % just put at end
            altFormula = sprintf('%s + %s', nullFormula, altFormula);
    
        else % if there are random effects, put new term before they start
    
            % put just after before the brackets
            altFormula = sprintf('%s %s + %s', nullFormula(1:parenInd-1), altFormula, nullFormula(parenInd:end));
        end
    end
    
    if ~exist('fitFunc','var') || isempty(fitFunc)
        % pick simplest one that will work
    
        % are there random effects?
        isRE = ~isempty(regexp(nullFormula, '(','ONCE'));
        isGLM = ~isempty(glmArgs); % are there glmeArguments
    
        if isRE
            if isGLM
                fitFunc = @fitglme;
            else
                fitFunc = @fitlme;
            end
        else
            if isGLM
                fitFunc = @fitglm;
            else
                fitFunc = @fitlm;
            end
        end
    end
    
    % refit() is quicker than fitting a new model, but only works for fitglme
    if strcmp(char(fitFunc),'fitglme')
        useRefit = 1; % ~25% faster to use this
    else
        useRefit = 0; % no refit() method for these functions
    end
    
    
    %% do the base regressions and check they differ by one term only
    
    nullFit = fitFunc(dataTab, nullFormula, glmArgs{:});
    altFit = fitFunc(dataTab, altFormula, glmArgs{:});
    
    altInNull = ismember(altFit.CoefficientNames, nullFit.CoefficientNames);
    
    if sum(~altInNull) ~= 1
        error('altFormula must contain ONE term that is not in the nullFormula');
    end
    
    altTermInd = find(~altInNull); % index of altTerm
    
    % print term using
    fprintf('running parametric bootstrap on the term %s, with %d bootstraps, and function %s\n', altFit.CoefficientNames{altTermInd}, nBoots, char(fitFunc));
    
    trueBeta = altFit.Coefficients.Estimate(altTermInd); % get true effect in data
    origP = altFit.Coefficients.pValue(altTermInd); % and the original p-value
    
    %% do bootstrapping
    
    
    regTab1 = dataTab; % copy this to overwrite data
    nullBetas = zeros(nBoots,1); % preallocate
    
    for i = 1:nBoots
        if mod(i,100)==0; disp(i); end % counter
    
        if useRefit % only works for glme, faster than other method
            altFit1 = refit(altFit, random(nullFit));
        
        else % for lme/glm/lm must git a new model
            
            % generate samples from the null model
            regTab1.(nullFit.ResponseName) = random(nullFit, dataTab);
    
            % fit altModel
            altFit1 = fitFunc(regTab1, altFormula, glmArgs{:});
        end

        
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


