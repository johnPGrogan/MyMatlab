function [corrP, t_obs, betas, se, df] = lmeEEG_categorical(eegMatrix, behTab, formula, nPerms, tail, chan_hood, categories, contrasts)
% function [corrP, t_obs, betas, se, df] = lmeEEG_categorical(eegMatrix, behTab, formula, nPerms, tail, chan_hood, categories, contrasts)
% 
% lmeEEG for categorical predictors - marginal EEG is done the same, but
% now we get a t-stat for each categorical level, and can also get a F-stat
% for the overall effect (stored as additional final row in corrP and
% t_obs). Since F-stats are always positive, and FindClusters only looks
% for negative t-stats (assumes symmetrical), we will invert all F-stats
% and only do a one-tailed comparison on the F-distribution. Can pass in
% cell array of contrasts to test also, which will be permuted with F-tests
% too
% 
% Call lmeEEG pipeline on one-channel's data - quicker way to do
% mixed-effects permutation testing, by first regressing out random
% effects.
% It regresses out random-effects at each time-point, and rebuilding the
% Fixed-effects + residuals as the 'marginal effects matrix', and then it
% runs a fixed-effects regression on those to get the "true" t-values. It
% then runs a permutation test, swapping conditions within each
% participant, and regressing the fixed-effects to build the null
% distribution, and uses my FindClusterLikeGND.m (based on DMGroppe toolbox)
% to get the clusters, rather than lmeEEG's TFCE as that does not seem to
% work for 1-channel.
% 
% Uses parfor, col, nanzscore, FindClusterLikeGND, and lmeEEG toolbox
% 
% 
% Inputs:
%   eegMatrix = EITHER:
%           1. matrix of size [nPP, nTimes, nLevels, nTr] (1 channel only), OR
%           2. matrix of size [nPP*nTr, nTimes, nChannels] - in which case,
%              you must pass in behTab too!!
%   behTab = [OPTIONAL] table with behavioural data and 'pp1' in, if passed
%       in, assumes that eegMatrix is a matrix or table of same height.
%       factor used in formula should be categorical with reference level
%       set
%   formula = [OPTIONAL] regression formula with Random Effects (DV must be 'EEG'),
%     and fixed-effects either just '1 + fac' if not using behTab, or
%     variables in behTab if that is given. default is 'EEG ~ 1 + fac + (1 | pp1)'
%   nPerms = number of permutations to run
%   tail = [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference). F-stat permutation for overall effect
%       across categorical levels must be one-tailed as F-stats are
%       positive, so that will be forced, but other predictors will use
%       tail.
%   chan_hood = if 1 channel given, is set to false, otherwise should be:
%     chan_hood = spatial_neighbors(eeg.chanlocs, 0.61, []);
%   categories = if not giving behTab, then categories can be given to set
%     the order of categories (1st one becomes reference level) to match
%     1:nL
%   contrasts = cell array of contrast matrices (bools), allows multiple 
%     rows, 1 column per, FE term in model. default is none, or can pass in
%     {'overall'} to just test that.
% 
% Outputs (no longer returns intercept row):
%   corrP = cluster-corrected p-values [nCoeff, nTimes]
%   t_obs = [nCoeff nTimes] "true" t-values from the 'marginal effects
%     matrix', i.e., with the RandomEffects regressed out
%   betas = [nCoeff, nTimes] "true" beta coeffs from marginal matrix
%   se = [nCoeff, nTimes] standard errors of "true" beta coeffs
%   df = degrees of freedom 
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

if ~exist('chan_hood','var') || isempty(chan_hood)
    chan_hood = false; % one channel
end

if ~exist('contrasts','var') || isempty(contrasts)
    contrasts = {}; % default is none, not even overall
    nC = 0; % 
else
    nC = length(contrasts);
    if nC==1 && ischar(contrasts{1}) % if 1, can be just {'overall'} 
        assert(strcmpi(contrasts{1}, 'overall'), "if contrasts is string, must be {'overall'}");
    else % otherwise must be [0 1 1] format
        assert(all(cellfun(@islogical, contrasts) | cellfun(@isnumeric, contrasts)),...
            "multiple contrasts must be cell array of {[0 1 0]; [0 0 1]} format for multiples")
    end

    fprintf('\nTesting following %d contrasts:', nC);
    contrasts{:}
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

    % make fac categorical
    warning('making fac column into a categorical, with 1st category as reference');
    dataTab.fac = categorical(dataTab.fac);

    nNonNan = sum(~isnan(eegMatrix(:,1,:,:)),'all'); % non-nan trials for 1st timepoint
end

% check sizes
assert(size(dataTab,1) == size(dvMat,1), 'dataTab and dimension mismatch');


%% need to make sure there are no NaNs

catCol = table2array(varfun(@iscategorical, dataTab));

if any(catCol)
    fprintf('\ncategorical predictor is: %s\n', dataTab.Properties.VariableNames{catCol});
else
    error('not categorical predictor');
end

% categorical can't be nan, or be table2array, skip for now
% remove all-nan rows?
toRemove = all(isnan(dvMat),[2 3]) | all(isnan(table2array(dataTab(:,~catCol))),2);

assert(sum(~toRemove) == nNonNan, 'NaN mismatch');

dataTab(toRemove,:) = [];
dvMat(toRemove,:,:) = [];


% if any(isnan(table2array(dataTab)),'all')
%     % run one regression to get predictor names, and then discard rows with
%     % nans in those
%     dataTab1 = dataTab;
%     dataTab1.EEG = dvMat(:,1,1);
%     m = fitglme(dataTab1, formula);
%     colNames = m.Formula.PredictorNames; % get names to keep
%     dataTab = dataTab(:, ismember(dataTab.Properties.VariableNames, ['EEG', colNames])); % remove other columns for parfor
% 
%     toRemove = any(isnan(table2array(dataTab)),2); % remove any nan rows from this
%     dataTab(toRemove,:) = [];
%     dvMat(toRemove,:,:) = [];
% 
%     clearvars dataTab1 m colNames;
% end

%% zscore predictors after removals

v = dataTab.Properties.VariableNames; % store
dataTab(:, ~catCol) = varfun(@nanzscore, dataTab(:,~catCol)); % zscore each column after removing NaNs
dataTab.Properties.VariableNames = v; % replace names

%% set categorical order, 1st is refence level

if exist('categories','var') && ~isempty(categories)
    fprintf('\n Re-ordering categories');
    dataTab.(v{catCol}) = reordercats(dataTab.(v{catCol}), categories);
end

%% 

[nRows, nT, nCh] = size(dvMat);
fprintf('\n%d Rows, %d Time-points, %d Channels, %d Permutations', nRows, nT, nCh, nPerms)

%% regress out RE, leaving just fitted FE + residuals

[mEEG, m1] = regressOutRE(dataTab, dvMat, formula);
X = designMatrix(m1);
nFE = size(X,2); % number of fixed effects + 1 for f-stat

df = [nFE-1,  m1.DFE]; % degrees of freedom for later [effect, resid]

clear dataTab; % reduce memory


%% get 'true' FE effects from this marginal data

fprintf('\nRunning regressions on marginals');
[t_obs, betas, se] = deal(NaN(nFE+nC, nT, nCh));
parfor iCh = 1:nCh
    EEG = mEEG(:,:,iCh); % copy
    [t_obs(:,:,iCh), betas(:,:,iCh), se(:,:,iCh)] = lmeEEG_regress_withF(EEG, X, contrasts);
end


%% permutation test

% this only returns unique ones, so if few trials per person it can get
% stuck
fprintf('\nCreating row permutations');
[rowPerms] = lmeEEG_permutations_repeats(m1.Variables.pp1, nPerms); % nperms within-subjects permutations of X
% this version allows non-unique permutations within a person, but will
% still given unique permutations across entire dataset if a lot of trials

fprintf('\nRunning %d permutations: ', nPerms);
t_perms = NaN(nFE+nC,nT,nCh,nPerms); % Initialize t-map
parfor iP = 1:nPerms
    XX = X(rowPerms(:,iP),:); % get indices for this perm
    if mod(iP,nPerms/10)==0; disp(iP); end % fprintf does not get output within parfor, only at end, so use disp
        
    for iCh = 1:nCh
%         EEG = mEEG(:,:,iCh); % copy   
        [t_perms(:,:,iCh,iP)] = lmeEEG_regress_withF(mEEG(:,:,iCh), XX, contrasts); % this works across row of samples now
    end
end

%% findclust
% for t-statistics

fprintf('\nFinding clusters: ');

n = nFE-1; % ignore intercept
corrP = NaN(n+nC, nT, nCh); % will skip intercept
for i = 1:n % skip intercept
    % make inputs [nT nCh (nPerms)]
    corrP(i,:,:) = FindClustersLikeGND(shiftdim(t_obs(i+1,:,:),1), shiftdim(t_perms(i+1,:,:,:),1), chan_hood, tail, df(2)); %[times chans]
end

%% find cluster for F-stats
% as F >=0, need to use different version of FindClusters for each contrast

for i = 1:nC
    corrP(i+n,:,:) = FindClustersLikeGND_F(shiftdim(t_obs(i+1+n,:,:),1), shiftdim(t_perms(i+1+n,:,:,:),1), chan_hood, 1, df); %[times chans]
end
%% remove intercepts from returned values

% t_perms(1,:,:,:) = [];
t_obs(1,:,:) = [];
betas(1,:,:) = [];
se(1,:,:) = [];
