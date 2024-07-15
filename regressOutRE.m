function [mEEG, m1] = regressOutRE(dataTab, dvMat, formula)
% function [mEEG, m1] = regressOutRE(dataTab, dvMat, formula)
% 
% To prepare to run lmeEEG_regress, first regress out random-effects from
% each channel + time-point, returning just the Fixed-effects + Residuals
% in a matrix the same size as dvMat.
% 
% Uses parfor across time-points
% 
% Inputs:
%   dataTab: [nPP*nTr, ...] table with columns for FE and RE
%   dvMat: [nPP*nTr, nTimes, nChans] matrix of EEG values - should not
%     contain NaN values
%   formula: regression formula, with DV name as 'EEG'. Must be the same
%     formula you will use for the lmeEEG regressions
% 
% Outputs:
%   mEEG = 'marginal effects matrix' of Fixed-effects + residuals (i.e.
%     random effects removed).
%   m1 = regression fitted to the 1st timepoint + channel, so has the
%     design matrix, coefficient names, and degrees of freedom
% 
% 


[~, nT, nCh] = size(dvMat);

if isa(dvMat, 'single'); dvMat = double(dvMat); end % fitlme needs double not single

%% run one regression first to get the FE names and design matrix

dataTab.EEG = double(squeeze(dvMat(:,1,1))); 
t1=tic;
m1 = fitlme(dataTab, formula);
t1 = toc(t1);
fprintf('\nFitting one time/channel took %g seconds', t1);
% X = designMatrix(m1);
% nFE = size(X,2); % number of fixed effects
% df = m1.DFE; % degrees of freedom for later

colNames = m1.Formula.PredictorNames; % get names to keep
eegTab = dataTab(:, ismember(dataTab.Properties.VariableNames, ['EEG', colNames])); % remove other columns for parfor

%% now regress each time+chan

fprintf('\nBuilding marginal effects, removing random-effects with formula:\n %s', formula );
mEEG = single(NaN(size(dvMat)));
n=nT*nCh; % parfor across both time*chans at once to be quicker
% for iCh = 1:nCh
%     if iCh>1; fprintf('%d, ', iCh); end
    parfor i = 1:n
        if mod(i, round(n/10))==0; disp(i); end
        eegTab1 = eegTab; % copy
        eegTab1.EEG = dvMat(: ,i); % overwrite EEG

        m = fitlme(eegTab1, formula); % fit it

        mEEG(:,i) = single(fitted(m,'Conditional',0) + residuals(m)); % Extract marginal EEG = FE + residuals
    end
% end