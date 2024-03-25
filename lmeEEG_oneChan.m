function [corrP, t_perms, t_obs, betas, se, df] = lmeEEG_oneChan(eegMatrix, behTab, formula, nPerms, tail, chan_hood)
% function [corrP, t_perms, t_obs, betas, se, df] = lmeEEG_oneChan(eegMatrix, behTab, formula, nPerms, tail, chan_hood)
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
%       in, assumes that eegMatrix is a matrix or table of same height
%   formula = [OPTIONAL] regression formula with Random Effects (DV must be 'EEG'),
%     and fixed-effects either just '1 + fac' if not using behTab, or
%     variables in behTab if that is given. default is 'EEG ~ 1 + fac + (1 | pp1)'
%   nPerms = number of permutations to run
%   tail = [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference)
%   chan_hood = if 1 channel given, is set to false, otherwise should be:
%     chan_hood = spatial_neighbors(eeg.chanlocs, 0.61, []);
% 
% Outputs:
%   corrP = cluster-corrected p-values [nCoeff, nTimes] (1st coeff is intercept)
%   t_perms = [nCoeff, nTimes, 1, nPerms] t-values for each coefficient at
%     each permutation
%   t_obs = [nCoeff nTimes] "true" t-values from the 'marginal effects
%     matrix', i.e., with the RandomEffects regressed out
%   betas = [nCoeff, nTimes] "true" beta coeffs from marginal matrix
%   se = [nCoeff, nTimes] standard errors of "true" beta coeffs
%   df = degrees of freedom 
%  

if ~exist('nPerms','var')
    nPerms = 1000;
end
if ~exist('tail','var')
    tail = 0; % two-tailed
end

if ~exist('formula', 'var')
    formula = 'EEG ~ 1 + fac + (1 | pp1)'; % RE will be removed
end

if ~exist('chan_hood','var')
    chan_hood = false; % one channel
end

%% can pass in eeg table? and beh table? and then FE + RE?
if exist('behTab','var')
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

v = dataTab.Properties.VariableNames; % store
dataTab = varfun(@nanzscore, dataTab); % zscore each column after removing NaNs
dataTab.Properties.VariableNames = v; % replace names

[~, nT, nCh] = size(dvMat);
fprintf('\n%d Time-points, and %d Channels', nT, nCh)

%% regress out RE, leaving just fitted FE + residuals

[mEEG, m1] = regressOutRE(dataTab, dvMat, formula);
X = designMatrix(m1);
nFE = size(X,2); % number of fixed effects
df = m1.DFE; % degrees of freedom for later

clear dataTab; % reduce memory


%% get 'true' FE effects from this marginal data

fprintf('\nRunning regressions on marginals');
[t_obs, betas, se] = deal(NaN(nFE, nT, nCh));
parfor iCh = 1:nCh
    EEG = mEEG(:,:,iCh); % copy
    [t_obs(:,:,iCh), betas(:,:,iCh), se(:,:,iCh)] = lmeEEG_regress(EEG, X)
end


%% permutation test

% this only returns unique ones, so if few trials per person it can get
% stuck
fprintf('\nCreating row permutations');
[rowPerms] = lmeEEG_permutations_repeats(m1.Variables.pp1, nPerms); % nperms within-subjects permutations of X
% this version allows non-unique permutations within a person, but will
% still given unique permutations across entire dataset if a lot of trials

fprintf('\nRunning %d permutations: ', nPerms);
t_perms = NaN(nFE,nT,nCh,nPerms); % Initialize t-map
parfor iP = 1:nPerms
    XX = X(rowPerms(:,iP),:); % get indices for this perm
    if mod(iP,nPerms/10)==0; disp(iP); end % fprintf does not get output within parfor, only at end, so use disp
        
    for iCh = 1:nCh
%         EEG = mEEG(:,:,iCh); % copy   
        [t_perms(:,:,iCh,iP)] = lmeEEG_regress(mEEG(:,:,iCh), XX); % this works across row of samples now
    end
end

%% findclust

fprintf('\nFinding clusters');

corrP = NaN(nFE, nT, nCh);
for i = 1:nFE
    % make inputs [nT nCh (nPerms)]
    corrP(i,:,:) = FindClustersLikeGND(shiftdim(t_obs(i,:,:),1), shiftdim(t_perms(i,:,:,:),1), chan_hood, tail, df); %[times chans]
end


