function [p, trueT, trueB, origP, nullT] = ParametricBootstrapTime(dataTab, responseVars, fullFormula, termOfInterest, nBoots, glmArgs, tail, doPlot, fitFunc)
% function [p, trueT, trueB,origP, nullT] = ParametricBootstrap(dataTab, fullFormula, termOfInterest, nBoots, glmArgs, tail, doPlot, fitFunc)
%
% Perform a parametric bootstrap (Bůžková, Lumley & Rice, 2011) in order to
% see how likely an interaction term-of-interest is under a null model. 
% This function applies this across a range of responses (e.g. times or
% channels), to keep the FWER across the whole time-series at .05. We 
% fit the null model (which is only missing the term-of-interest) to the
% data, and then generate samples from this null model many times,
% fitting the alternative model and keeping the coefficient estimate for
% the interaction term. This gives a distribution of the interaction from
% the null model which does not contain it.
% We do this for each responseVar (timepoint).
% We then calculate the p-value as the probability of the true estimate
% given this null distribution.
%
%
% 
% Inputs:
%   dataTab = table to pass to fitglme for regression, each term in the
%       formulae must be a column. Variables should be z-scored
%   responseVars = cell array of response variable names [1 nRespVars]
%   fullFormula = the full formula with all terms (DV should be '%s'), 
%       including the term-of-interest and any random effects. The 
%       term-of-interest will be removed from this to create the null model.
%       e.g. '%s ~ 1 + x1 + x2 + x1:x2 + (1 | pp)'
%   termOfInterest = one of the terms from the fullFormula which will be
%       removed to create the null model, and sampled from this. If an
%       interaction, it should be in the form 'a:b' not 'a*b', e.g. if the
%       fullFormula = '%s ~ 1 + a*b', this will be expanded into '%s ~ 1 + a
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
%   p = p-value (probability of the alternative term-of-interest under the
%       null distribution
%   trueT = alternate term-of-interest t-statistic from the
%     alterate model fit to the real data (i.e. the 'real' t-stat) [1 nRespVars]
%   trueB = alternate term-of-interest coefficient estimate from the
%     alterate model fit to the real data (i.e. the 'real' estimate) [1 nRespVars]
%   origP = alternate term-of-interest p-value from fitting the alternate
%       model to the read data (i.e. 'real' p-value) [1 nRespVars]
%   nullT = coefficient estimates of the alternate model
%     term-of-interest, generated from the null model [nBoots x nRespVars]
%
% John Grogan, 2021



    %% check inputs
    
    if ~istable(dataTab); error('dataTab must be a table'); end

    if ~iscell(responseVars)
        if ischar(responseVars)
            responseVars = {responseVars}; % may just be one string for one var
        else
            error('responesVars should be a cell array of strings');
        end
    else
        nRV = length(responseVars);
    end

    
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
    
    iRV = 1; % respo
    altFit = fitFunc(dataTab, sprintf(fullFormula, responseVars{iRV}), glmArgs{:}); % use 1st response var
        
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

        nullFormula = ['%s ~ ', terms]; % keep %s at beginnig
    end


    % get the null model
    nullFit = fitFunc(dataTab, sprintf(nullFormula, responseVars{iRV}), glmArgs{:});
    
    % check they only differ by 1 FE term
    altInNull = ismember(altFit.CoefficientNames, nullFit.CoefficientNames);
    
    if sum(altInNull==0) ~= 1
        error('altFormula must contain ONE term that is not in the nullFormula:\n nullFormula: %s \nfullFormula: %s\n', nullFormula, fullFormula);
    end
    
    altTermInd = find(~altInNull); % index of altTerm
    
    % print term using
    fprintf('\nrunning parametric bootstrap on the term: %s, \nnullModel: %s \nFullModel: %s \nNum bootstraps: %d (maximum p-value resolution = %g), \nFitting function: %s \nTail: %d\n', ...
        termOfInterest, nullFormula, fullFormula, nBoots, 1/nBoots, char(fitFunc),tail);

    
    %% do bootstrapping
    
    
%     regTab1 = dataTab; % copy this to overwrite data
    [nullT, nullB] = deal(zeros(nBoots,nRV)); % preallocate
    disp('starting bootstraps');
    parfor iRV = 1:nRV
        fprintf('\niRV=%d: i=', iRV);
        regTab1 = dataTab;
        altFit = fitFunc(dataTab, sprintf(fullFormula, responseVars{iRV}), glmArgs{:}); % fit full model, so can use refit later
        trueB(1,iRV) = altFit.Coefficients.Estimate(altTermInd); % get true effect in data
        origP(1,iRV) = altFit.Coefficients.pValue(altTermInd); % and the original p-value
        trueT(1,iRV) = altFit.Coefficients.tStat(altTermInd); % get true effect in data

        nullFit = fitFunc(dataTab, sprintf(nullFormula, responseVars{iRV}), glmArgs{:}); % get null model to simulate data from

        for i = 1:nBoots
            if mod(i,100)==0; fprintf('%d, ',  i); end % counter
        
            if useRefit % only works for glme, faster than other method
                altFit1 = refit(altFit, random(nullFit));
            
            else % for lme/glm/lm must git a new model
                
                % generate samples from the null model
                regTab1.(nullFit.ResponseName) = random(nullFit, dataTab);
        
                % fit altModel
                altFit1 = fitFunc(regTab1, sprintf(fullFormula, responseVars{iRV}), glmArgs{:});
            end
    
            
            % get beta
            nullB(i,iRV) = altFit1.Coefficients.Estimate(altTermInd);
            nullT(i,iRV) = altFit1.Coefficients.tStat(altTermInd);
%             a = anova(altFit1);
%             f(i,1) = a.FStat(altTermInd);
        end
    end
    

    %% take only most extreme per time?

%     [m,i] = max(abs(nullT),[],2);
    
%     nullT = nullT(:);% use all samples generated across all times

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
%         p = 1 - abs(mean(nullT <= trueT) - mean(nullT >= trueT)); % each time sep
        p = 1 - abs(mean(nullT(:) <= trueT) - mean(nullT(:) >= trueT)); % all time samples together
%         p = 1 - abs(mean(min(nullT,[],2) <= trueT) - mean(max(nullT,[],2) >= trueT)); % using min/max method
    end
    
    %% estimate given FWER
    
    if tail == 0 % two-tailed
%         bMax = prctile(abs(nullT(:)), 95);
%         alphaEst = mean(abs(nullT(:)) >= bMax);
        bMax = prctile(nullT(:), [2.5 97.5]);
        alphaEst = mean(~isBetween(nullT(:), bMax,0));

    elseif tail == 1 % one tailed
        bMax = prctile(nullT(:), 97.5);
        alphaEst = mean(nullT(:) >= bMax);
    elseif tail == -1 % one tailed
        bMax = prctile(nullT(:), 2.5);
        alphaEst = mean(nullT(:) <= bMax);
    end
    fprintf('\nEstimated FWER: %f\n',alphaEst);
end


