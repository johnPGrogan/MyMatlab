function statsTable = myStatsExtract(glmeFits, colNames, contrasts)
% function myStatsExtract(glmeFits, colNames, contrasts)
% will extract the p-values and betas from each glmeFit, and also do an anova
% and get the p&F for the main effects that are not in glmeFit
% and if contrasts is given, will get the p&F for each contrast
% works on a cell array of glmes 
% 
% Inputs:
%   glmeFits = cell array of fitglme outputs, each with same effects but
%     different DVs
%   colNames = cell array of names of DVs in glmeFits
%   contrassts = contrasts to test, corresponding to coefficients in
%     glmeFits
% 
% Outputs:
%  statsTable = table with 1 column per DV. rows are [p; F ] for the
%    anova main effects (for those which do not appear in glmeFits - i.e.
%    only to test overall effects of categorical contrasts), followed by
%    the main effects of coefficients from the lme [p; Beta], then the
%    contrats [p;F]. Intercepts are not extracted
% 

if ~exist('contrasts','var') || isempty(contrasts)
    nC = 0;
else
    nC = size(contrasts,1);
end

statsTable = table; % set rows

nF = length(glmeFits);

for i = 1:nF
    % do anova on this to get the overall effect across levels
    a = anova(glmeFits{i});
    % remove main effects that are the same as in lme, including intercept
    a(ismember(a.Term, glmeFits{i}.CoefficientNames'),:) = [];



    coefStats = NaN(nC,2);
    for j = 1:nC
        [coefStats(j,1), coefStats(j,2)] = coefTest(glmeFits{i},contrasts(j,:)); % get pos vs neg
    end


    % combine it in order [p; F or Beta] for each test
    statsTable.(colNames{i}) = col([a.pValue, a.FStat;... % anova of main effects
        glmeFits{i}.Coefficients.pValue(2:end), glmeFits{i}.Coefficients.Estimate(2:end);... % lme for each contrast
        coefStats;... % comparison of other contrasts
        ]');

end

% get row-names from the stats
r = [a.Term; glmeFits{i}.CoefficientNames(2:end)']; % get anova + lme names

% % for contrasts, extract the names and store difference
% b = glmeFit.CoefficientNames(contrasts~=0)'; % get contrasts - difference as name
% len = cellfun(@length, b); [~,len] = sort(len); % which is shorter
% r{end+1} = '';
% for j = 1:length(b{len(1)})
%     if b{1}(j)==b{2}(j)
%         r{end}(end+1) = b{1}(j);
%     else
%         break
%     end
% end
% r{end} = [r{end} b{1}(j:end) '_vs' b{2}(j-1:end)];
% 
% % would only work for 2 contrasts
% % instead, just store the contrast as a string
for j=1:nC
    r{end+1} = sprintf('%d_',contrasts(j,:));
end

r2 = [repmat({'F_'}, size(a,1),1); repmat({'B_'}, glmeFits{i}.NumCoefficients-1,1); repmat({'F_'}, nC, 1)];
statsTable.Properties.RowNames = strcat(col([repmat({'p_'},size(r,1),1), r2]'), col(repmat(r,1,2)'));
