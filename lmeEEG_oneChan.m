function [corrP, t_perms, t_obs, betas, se, df] = lmeEEG_oneChan(eegMatrix, nPerms, tail)
% function [corrP, t_perms, t_obs, betas, se, df] = lmeEEG_oneChan(eegMatrix, nPerms, tail)
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
%   eegMatrix = [nPP nTimes nLevels nTr] (is less sensitive if using
%     aggregate-level (i.e. averaged over trials)
%   nPerms = number of permutations to run
%   tail = [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference)
% 
% Outputs:
%   corrP = cluster-corrected p-values [nCoeff, nTimes] (1st coeff is intercept)
%   t_perms = [nCoeff, nTimes, 1, nPerms] t-values for each coefficient at
%     each permutation
%   t_obs = [nCoeff nTimes] "true" t-values from the 
%  

if ~exist('nPerms','var')
    nPerms = 1000;
end
if ~exist('tail','var')
    tail = 0; % two-tailed
end

chan_hood = false; % one channel

% can pass in eeg table? and beh table? and then FE + RE?

%% set up matrix + table

[nPP, nT, nL, nTr] = size(eegMatrix);

dvMat = reshape(permute(eegMatrix,[1 3 4 2]),nPP*nL*nTr, nT); %[nPP*nL*nTr, nT]

dataTab = [col(repmat(1:nPP, 1, nL*nTr)), col(repmat(1:nL, nPP,nTr)), col(repmat(1:nTr,nL*nPP,1))]; %[pp factor tr]

% dataTab = nanzscore(dataTab); % mean-centre and standardise
% better to do this after removing NaN rows?

dataTab = array2table(dataTab, 'VariableNames', {'pp1','fac','trial'}); % trial is unused


%% need to make sure there are no NaNs

% remove all-nan rows?
toRemove = all(isnan(dvMat),2);

assert(sum(~toRemove) == sum(~isnan(eegMatrix(:,1,:,:)),'all'), 'nan mismatch');

dataTab(toRemove,:) = [];
dvMat(toRemove,:,:) = [];

v = dataTab.Properties.VariableNames; % store
dataTab = varfun(@nanzscore, dataTab); % zscore each column
dataTab.Properties.VariableNames = v; % replace names

[~, nT, nCh] = size(dvMat);

%% regress out RE, leaving just fitted FE + residuals

formula = 'EEG ~ 1 + fac + (1 | pp1)'; % RE will be removed

pp1 = nominal(dataTab.pp1); % same as categorical
fac = dataTab.fac; % in tutorial, this was categorical too, though they had only 2 levels,

fprintf('\nBuilding marginal effects, removing random-effects with formula:\n %s', formula );
mEEG = nan(size(dvMat));
for iCh = 1:nCh
    parfor iT = 1:nT
        EEG = dvMat(:,iT,iCh);
        EEG = table(EEG, fac, pp1);

        m = fitlme(EEG, formula); % fit it

        mEEG(:,iT,iCh) = fitted(m,'Conditional',0) + residuals(m); % Extract marginal EEG = FE + residuals
    end
end

% Extract design matrix X
EEG = double(squeeze(dvMat(:,1,1)));
EEG = table(EEG, fac, pp1);
m = fitlme(EEG, formula);

X = designMatrix(m);
nFE = size(X,2);

%% get 'true' FE effects from this marginal data

fprintf('\nRunning regressions on marginals');
[t_obs, betas, se] = deal(NaN(nFE, nT, nCh));
for iCh = 1:nCh
    parfor iT = 1:nT
        EEG = mEEG(:,iT,iCh);
        [t_obs(:,iT,iCh), betas(:,iT,iCh), se(:,iT,iCh)] = lmeEEG_regress(EEG, X)
    end
end


%% permutation test

% this only returns unique ones, so if doing aggregate, it gets stuck
[rowPerms] = lmeEEG_permutations(pp1,nPerms); % nperms within-subjects permutations of X

fprintf('\nRunning %d permutations:', nPerms);
t_perms = NaN(nFE,nT,nCh,nPerms); % Initialize t-map
for iP = 1:nPerms
    if mod(iP,50)==0; fprintf('%d, ', iP); end

    XX = X(rowPerms(:,iP),:); % get indices for this perm
    for iCh = 1:nCh
        parfor iT = 1:nT
            EEG = squeeze(mEEG(:,iT,iCh));
            [t_perms(:,iT,iCh,iP)] = lmeEEG_regress(EEG, XX);
        end
    end
end


%% findclust

fprintf('\nFinding clusters');
df = m.DFE; % degrees of freedom

corrP = NaN(nFE, nT, nCh);
for i = 1:nFE
    % make inputs [nT nCh (nPerms)]
    corrP(i,:,:) = FindClustersLikeGND(shiftdim(t_obs(i,:,:),1), shiftdim(t_perms(i,:,:,:),1), chan_hood, tail, df)'; %[times chans]
end


