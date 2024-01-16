function coeffs = standardiseRegCoeffs(fit)
% returns standardised regression coefficients, SE, upper + lower CI
% does this by applying conversion coefficient std(Xi)/std(Y) for
% fitglme model with regression Y ~ 1 + X(1) + ... + X(n)
% It does not scale the intercept, and will not work if there are any 
% interactions! as these change the relationships between the variables, 
% they are not simply scaled.
%
% Inputs:
%   fit = model fit from fitlm, fitglm, fitlme, fitglme. Contains 
%       'Coefficients' field which is a table (or lme titleddataset)
% 
% Outputs:
%  coeffs = 'Coefficients' field which is the table, with the following
%       columns standardised: Estimate, SE, Upper, Lower (i.e. upper/lower
%       Confidence Intervals). (Upper + Lower only present in mixed-effects
%       model fits using fitlme or fitglme). t-values, p-values, and DF are
%       unchanged.


%% get coeffs

feTerms = fit.CoefficientNames;
feTerms(ismember(fit.CoefficientNames,'(Intercept)')) = []; % remove intercept FE

%% check no interactions

isInteraction = any(cellfun(@(x) sum(ismember(x,':')), feTerms)); % max level in FE-terms

assert(~isInteraction, 'interaction term present, cannot simply convert. instead standardise your variables');

nFE = length(feTerms);
assert(nFE, 'no fixed effects in model');

%% get current coeffs
coeffs = fit.Coefficients; 
coeffs(ismember(fit.CoefficientNames,'(Intercept)'),:) = []; % remove intercept


if regexp(class(fit), 'MixedModel')
    
    %% are there random effects?
    % if other than an intercept, give a warning that these will not be
    % standardised, and there could be differences between this standardising
    % and regressing standardised variables if the standardising changes the
    % relationships
    
    nRE = length(fit.Formula.RELinearFormula); % number of brackets
    nREs = zeros(nRE,1);
    for i = 1:nRE
        nREs(i) = fit.Formula.RELinearFormula{i}.NTerms;
        
        if any(~ismember(fit.Formula.RELinearFormula{i}.TermNames, '(Intercept)'))
            warning('random slopes found: %s, not converted as variable standardising changes these', fit.Formula.RELinearFormula{i}.LinearPredictor);
            warning('assumes the same relationship between variables (incl RE) with/out standardising');
        end
    end

    colNames = fit.Coefficients.Properties.VarNames; % get these

else
    colNames = fit.Coefficients.Properties.VariableNames;
end

%% create conversion factors

stdY = nanstd(fit.Variables.(fit.ResponseName)); % std of Y

stdX = NaN(nFE,1);
for i = 1:nFE
    stdX(i,1) = nanstd(fit.Variables.(feTerms{i}));
end

conversion = stdX ./ stdY;

%% apply

colsToConvert = ismember(colNames, {'Estimate', 'SE', 'Lower', 'Upper'});

for i = find(colsToConvert)
    coeffs.(colNames{i}) = coeffs.(colNames{i}) .* conversion;
end


