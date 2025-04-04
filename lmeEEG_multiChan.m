function [corrP, t_obs, betas, se, df, t_perms] = lmeEEG_multiChan(eegMatrix, behTab, formula, nPerms, tail, chanlocs, skipCluster, keepIntercept)
% function [corrP, t_obs, betas, se, df, t_perms] = lmeEEG_multiChan(eegMatrix, behTab, formula, nPerms, tail, chanlocs, skipCluster, keepIntercept)
% 
% Call lmeEEG pipeline for multi-channel data, uses TFCE instead of
% findCluster (quicker) and can skip that too if desired.
% It regresses out random-effects at each time-point, and rebuilding the
% Fixed-effects + residuals as the 'marginal effects matrix', and then it
% runs a fixed-effects regression on those to get the "true" t-values. It
% then runs a permutation test, swapping conditions within each
% participant, and regressing the fixed-effects to build the null
% distribution, and uses my lmeEEG_TFCE for clustering (takes a bit longer
% than FindClusers, but seems a bit more accurate?
% 
% Uses parfor, col, nanzscore, lmeEEG_TFCE, and lmeEEG toolbox
% 
% 
% Inputs:
%   eegMatrix = EITHER:
%           1. matrix of size [nPP, nTimes, nLevels, nTr] (1 channel only), OR
%           2. matrix of size [nPP*nTr, nTimes, nChannels] - in which case,
%              you must pass in behTab too!!
%   behTab = [OPTIONAL] table with behavioural data and 'pp1' in, if passed
%       in, assumes that eegMatrix is a matrix or table of same height
%   formula = [OPTIONAL] regression formula with Random Effects (DV must be 'EEG'),
%     and fixed-effects either just '1 + fac' if not using behTab, or
%     variables in behTab if that is given. default is 'EEG ~ 1 + fac + (1 | pp1)'
%   nPerms = number of permutations to run
%   tail = [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference)
%   chanlocs = eeg chanlocs structure
%   skipCluster = 1=don't find clusters, just return t_perms etc, default=0
%   keepIntercept = (default=0), 1=keep intercept and run through
%      clustering
% 
% Outputs (no longer returns intercept row):
%   corrP = cluster-corrected p-values [nCoeff, nTimes]
%   t_obs = [nCoeff nTimes] "true" t-values from the 'marginal effects
%     matrix', i.e., with the RandomEffects regressed out
%   betas = [nCoeff, nTimes] "true" beta coeffs from marginal matrix
%   se = [nCoeff, nTimes] standard errors of "true" beta coeffs from full
%     data - not from marginal data as those are too small
%   df = degrees of freedom 
%   t_perms = [nCoeffs nTimes nChans nPerms] permuted t-values
%

if ~exist('nPerms','var') || isempty(nPerms)
    nPerms = 1000;
end
if ~exist('tail','var') || isempty(tail)
    tail = 0; % two-tailed
end

if ~exist('formula', 'var') || isempty(formula)
    formula = 'EEG ~ 1 + fac + (1 + fac| pp1)'; % RE will be removed
end

if ~exist('skipCluster','var') || isempty(skipCluster)
    skipCluster = 0;
end
if ~skipCluster
    if ~exist('chanlocs','var') || isempty(chanlocs) % only needed if used
        error('chanlocs must be given if not skipping tfce clustering');
    end
    if ~exist('ept_ChN2','file')
        addpath(genpath('C:\Users\groganj1\OneDrive - TCDUD.onmicrosoft.com\GeneralScriptsOD\MatlabPackages\ept_TFCE-matlab-master'))
    end
end 
if ~exist('keepIntercept','var') || isempty(keepIntercept)
    keepIntercept = 0;
end
    
%% can pass in eeg table? and beh table? and then FE + RE?
if exist('behTab','var') && ~isempty(behTab)
    % if this is passed in, eegMatrix should be [nPP*nL*nTr, nT, nCh] matrix
    % and behTab should be [nPP*nTr, nBehVars] table

    dataTab = behTab; % copy [nPP*nTr, nbehvars]
    dvMat = eegMatrix; % copy [nPP*nTr, nTimes, nChans]

    nNonNan = sum(~isnan(eegMatrix(:,1,1,1))); % non-nan trials for 1st timepoint

else % no behTab, so eegMatrix is [nPP, nT, nL, nTr] for 1 channel
    %% set up matrix + table
    
    [nPP, nT, nL, nTr] = size(eegMatrix); % nTr is per level here
    
    dvMat = reshape(permute(eegMatrix,[1 3 4 2]),nPP*nL*nTr, nT); %[nPP*nL*nTr, nT]
    
    dataTab = [col(repmat(1:nPP, 1, nL*nTr)), col(repmat(1:nL, nPP,nTr)), col(repmat(1:nTr,nL*nPP,1))]; %[pp factor tr]
    
    % dataTab = nanzscore(dataTab); % mean-centre and standardise
    % better to do this after removing NaN rows?
    
    dataTab = array2table(dataTab, 'VariableNames', {'pp1','fac','trial'}); % trial is unused

    nNonNan = sum(~isnan(eegMatrix(:,1,:,:)),'all'); % non-nan trials for 1st timepoint
end

% check sizes
assert(size(dataTab,1) == size(dvMat,1), 'dataTab and dimension mismatch');


%% need to make sure there are no NaNs

% remove all-nan rows?
toRemove = all(isnan(dvMat),[2 3]) | all(isnan(table2array(dataTab)),2);

assert(sum(~toRemove) == nNonNan, 'NaN mismatch');

dataTab(toRemove,:) = [];
dvMat(toRemove,:,:) = [];

if any(isnan(table2array(dataTab)),'all')
    % run one regression to get predictor names, and then discard rows with
    % nans in those
    dataTab1 = dataTab;
    dataTab1.EEG = dvMat(:,1,1);
    m = fitglme(dataTab1, formula);
    colNames = m.Formula.PredictorNames; % get names to keep
    dataTab = dataTab(:, ismember(dataTab.Properties.VariableNames,  unique([m.Formula.ResponseName 'EEG', colNames]))); % remove other columns for parfor

    toRemove = any(isnan(table2array(dataTab)),2); % remove any nan rows from this
    dataTab(toRemove,:) = [];
    dvMat(toRemove,:,:) = [];

    clearvars dataTab1 m colNames;
end

%% zscore predictors after removals

v = dataTab.Properties.VariableNames; % store
dataTab = varfun(@nanzscore, dataTab); % zscore each column after removing NaNs
dataTab.Properties.VariableNames = v; % replace names

[nRows, nT, nCh] = size(dvMat);
fprintf('\n%d Rows, %d Time-points, %d Channels, %d Permutations', nRows, nT, nCh, nPerms)

%% regress out RE, leaving just fitted FE + residuals

[mEEG, m1, se] = regressOutRE(dataTab, dvMat, formula);
X = designMatrix(m1);
nFE = size(X,2); % number of fixed effects
df = m1.DFE; % degrees of freedom for later

clear dataTab; % reduce memory
clear dvMat;

%% get 'true' FE effects from this marginal data

fprintf('\nRunning regressions on marginals');
[t_obs, betas] = deal(single(NaN(nFE, nT, nCh)));
parfor iCh = 1:nCh
    EEG = mEEG(:,:,iCh); % copy
    [t_obs(:,:,iCh), betas(:,:,iCh)] = lmeEEG_regress(EEG, X)
end


%% permutation test

% this only returns unique ones, so if few trials per person it can get
% stuck
fprintf('\nCreating row permutations');
[rowPerms] = lmeEEG_permutations_repeats(m1.Variables.pp1, nPerms); % nperms within-subjects permutations of X
% this version allows non-unique permutations within a person, but will
% still given unique permutations across entire dataset if a lot of trials

fprintf('\nRunning %d permutations: ', nPerms);
t_perms = single(NaN(nFE,nT,nCh,nPerms)); % Initialize t-map
parfor iP = 1:nPerms
    XX = X(rowPerms(:,iP),:); % get indices for this perm
    if mod(iP,nPerms/10)==0; disp(iP); end % fprintf does not get output within parfor, only at end, so use disp
        
    for iCh = 1:nCh
%         EEG = mEEG(:,:,iCh); % copy   
        [t_perms(:,:,iCh,iP)] = lmeEEG_regress(mEEG(:,:,iCh), XX); % this works across row of samples now
    end
end


%% remove intercept, unless only one given

if m1.NumCoefficients > 1 && ~keepIntercept
    isInt = strcmp(m1.CoefficientNames, '(Intercept)'); % remove later, allows ' EEG ~ -1 + fac' to be run
    if any(isInt)
        fprintf('\nRemoving intercept from cluster-finding and outputs');
        t_perms(isInt,:,:,:) = [];
        t_obs(isInt,:,:) = [];
        betas(isInt,:,:) = [];
        se(isInt,:,:) = [];

        nFE = size(t_obs,1); % update
    end
end

%% findclust

if skipCluster
    fprintf('\nSkipping clusters:');
    corrP = single(NaN(size(t_obs)));
    return;
end

fprintf('\nFinding clusters:');

corrP = single(NaN(nFE, nT, nCh));
for i = 1:nFE % skip intercept
    % this uses parfor now
    tfceRes = lmeEEG_TFCE(permute(double(t_obs(i,:,:)),[3,2,1]), permute(double(t_perms(i,:,:,:)),[4,3,2,1]), chanlocs, [0.66 2]);
    corrP(i,:,:,1) = tfceRes.P_Values';

end


