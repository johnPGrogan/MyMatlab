function [p, trueT, trueB, origP, nullT] = ParametricBootstrap(dataTab, fullFormula, termOfInterest, nBoots, glmArgs, tail, doPlot, fitFunc)
% function [p, trueT, trueBeta, origP, nullBetas] = ParametricBootstrap(dataTab, fullFormula, termOfInterest, nBoots, glmArgs, tail, doPlot, fitFunc)
%
% Perform a parametric bootstrap (Bůžková, Lumley & Rice, 2011) in order to
% see how likely an interaction term-of-interest is under a null model. We 
% fit the null model (which is only missing the term-of-interest) to the
% data, and then generate samples from this null model many times,
% fitting the alternative model and keeping the coefficient estimate for
% the interaction term. This gives a distribution of the interaction from
% the null model which does not contain it.
% We then calculate the p-value as the probability of the true estimate
% given this null distribution.
%
%
% 
% Inputs:
%   dataTab = table to pass to fitglme for regression, each term in the
%       formulae must be a column. Variables should be z-scored
%   fullFormula = the full formula with all terms, including the
%       term-of-interest and any random effects. The term-of-interest will
%       be removed from this to create the null model.
%   termOfInterest = one of the terms from the fullFormula which will be
%       removed to create the null model, and sampled from this. If an
%       interaction, it should be in the form 'a:b' not 'a*b', e.g. if the
%       fullFormula = 'Y ~ 1 + a*b', this will be expanded into 'Y ~ 1 + a
%       + b + a:b', so termOfInterest should be 'a:b'.
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
%   p = p-value (probability of the alternate term-of-interest under the
%       null distribution)
%   trueT = t-statistic for the term-of-interest in the full model fit to
%       real data (i.e. the 'real' estimate)
%   trueB = alternate term-of-interest coefficient estimate from the
%     alterate model fit to the real data (i.e. the 'real' estimate)
%   origP = alternate term-of-interest p-value from fitting the alternate
%       model to the read data (i.e. 'real' p-value)
%   nullT = coefficient t-statistics for the alternate model
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

    
    if isempty(regexp(fullFormula, '~','ONCE')) % if it is not a full formula, just one term
        error('fullFormula should be a formula to pass into the lm');
    end


    if ~isempty(regexp(termOfInterest, '~','ONCE')) || ~isempty(regexp(termOfInterest, ' ','ONCE'))
        % if it is not a full formula, just one term
        error('termOfInterest should be one term with no spaces, not a formula');
    end
    
    if ~exist('fitFunc','var') || isempty(fitFunc)
        % pick simplest one that will work
    
        % are there random effects?
        isRE = ~isempty(regexp(fullFormula, '(','ONCE'));
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
    
    altFit = fitFunc(dataTab, fullFormula, glmArgs{:});
        
    if isempty(regexp(fullFormula,'\*','ONCE')) && ~isempty(regexp(fullFormula,termOfInterest,'ONCE')) % if not using ':', expand and add in RE
        % just pull term out of the given formula
        tInd = regexpi(fullFormula, termOfInterest); % start of it
        
        % need to also remove a plus either before or after
        plusInds = regexp(fullFormula,'\+');

        % which of those is closest to start/end of term?
        [~,i] = min(min(abs(plusInds - [tInd; tInd + length(termOfInterest)-1])));

        nullFormula = fullFormula; % copy
        nullFormula([plusInds(i) tInd:(tInd+length(termOfInterest)-1)]) = [];

        % check there aren't two consecutive +
        nullNoSpaces = regexprep(nullFormula,' ', '');
        if any(diff(regexp(nullNoSpaces,'+'))==1)
            error('nullFormula has too many spaces after removing the term of interest: %s', nullFormula);
        end
        

    else
        % remove the term of interest from the expanded model
        terms = altFit.CoefficientNames'; % get all fixed effect terms

        % replace 'intercept' with '1'
        interceptInd = strcmp(terms, '(Intercept)');
        if sum(interceptInd)==1
            terms(interceptInd) = {'1'};
        else
            terms = [{'-1'}; terms]; % if no intercept, must specify -1
        end    
       
        terms(strcmpi(terms, termOfInterest)) = []; % remove
        terms(:,2) = {' + '}; % to combine them
        terms = col(terms'); % change order
        terms = terms(1:end-1); % remove trailing plus
        terms = [terms{:}]; % make into one string

        reInd = regexp(fullFormula, '(','ONCE'); % are there random effects?
        if reInd
            terms = sprintf('%s + %s', terms, fullFormula(reInd:end)); % add them back in
        end

        nullFormula = sprintf('%s ~ %s', altFit.ResponseName, terms);
    end


    % get the null model
    nullFit = fitFunc(dataTab, nullFormula, glmArgs{:});
    
    % check they only differ by 1 FE term
    altInNull = ismember(altFit.CoefficientNames, nullFit.CoefficientNames);
    
    if sum(altInNull==0) ~= 1
        error('altFormula must contain ONE term that is not in the nullFormula:\n nullFormula: %s \nfullFormula: %s\n', nullFormula, fullFormula);
    end
    
    altTermInd = find(~altInNull); % index of altTerm
    
    % print term using
    fprintf('\nrunning parametric bootstrap on the term: %s, \nnullModel: %s \nFullModel: %s \nNum bootstraps: %d (maximum p-value resolution = %g) \nFitting function: %s \nTail: %d\n', ...
        termOfInterest, nullFormula, fullFormula, nBoots, 1/nBoots, char(fitFunc),tail);

    trueT = altFit.Coefficients.tStat(altTermInd); % get true t-stat in data
    origP = altFit.Coefficients.pValue(altTermInd); % and the original p-value
    trueB = altFit.Coefficients.Estimate(altTermInd); % and beta
    
    %% do bootstrapping
    
    
    regTab1 = dataTab; % copy this to overwrite data
    nullT = zeros(nBoots,1); % preallocate
    
    for i = 1:nBoots
        if mod(i,100)==0; disp(i); end % counter
    
        if useRefit % only works for glme, faster than other method
            altFit1 = refit(altFit, random(nullFit));
        
        else % for lme/glm/lm must git a new model
            
            % generate samples from the null model
            regTab1.(nullFit.ResponseName) = random(nullFit, dataTab);
    
            % fit altModel
            altFit1 = fitFunc(regTab1, fullFormula, glmArgs{:});
        end

        
        % get beta
        nullT(i) = altFit1.Coefficients.tStat(altTermInd);
        nullB(i) = altFit1.Coefficients.Estimate(altTermInd);
    end
    
    
    %% do plot?
    
    if doPlot
        figure();
        histogram(nullT);
        xline(trueT);
        xlabel('generated \beta from null model');
        ylabel('frequency');
    end
    
    %% calculate p-value
    
    
    if tail == 1
        p = mean(nullT >= trueT); % prob of getting as big or bigger
    
    elseif tail == -1
    
        p = mean(nullT <= trueT); % prob of getting as small or smaller
    
    elseif tail == 0
    
%         p = mean(abs(nullT) >= abs(trueT)); % prob of getting as extreme or more (i.e. pos and neg together)
        p = 1 - abs(mean(nullT <= trueT) - mean(nullT>= trueT));
    
    end


    %% estimate given FWER
    
    if tail == 0 % two-tailed
%         tMax = prctile(abs(nullT), 95);
%         alphaEst = mean(nullT >= tMax);
        tMax = prctile(nullT, [2.5 97.5]);
        alphaEst = mean(~isBetween(nullT, tMax,0));

    elseif tail == 1 % one tailed
        tMax = prctile(nullT, 97.5);
        alphaEst = mean(nullT >= tMax);
    elseif tail == -1 % one tailed
        tMax = prctile(nullT, 2.5);
        alphaEst = mean(nullT <= tMax);
    end
    fprintf('Estimated FWER: %f\n',alphaEst);

end


