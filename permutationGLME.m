function [corrP, trueT, trueB, trueP, permP, permT, permB] = permutationGLME(dvMat, dataTab, formula, coeffToPermute, tail, glmeArgs, chan_hood, nPerms, useCluster)
% function [corrP, trueT, trueB, trueP, permP, permT, permB] = permutationGLME(dvMat, dataTab, formula, coeffToPermute, tail, glmeArgs, chan_hood, nPerms, useCluster)
% Apply
% will be applied to every sample across time (or channel etc)
% and the p-values adjsuted via clustering to control FWER at .05 using
% FindClustersLikeGND() which uses the Mass Univariate Toolbox (DMGroppe),
% so that must be on the path.
%
% Inputs:
%   dvMat: matrix of data to use for DV, must match rows of dataTab [N time-samples chans]
%       should be zscored already if desired. fitglme will be run on each
%       column*page. If including channels, must supply chanLocs
%   dataTab: table with IVs and covars in, these can be categorical or
%       continuous
%   formula: regression formula for fitglme to use, with 'DV ~ ' at start
%   coeffToPermute: cell array of coefficient names to permute (must match
%     those in formula) - permutation is not properly defined for
%     interactions, so suggested to avoid those (and random-effects) or
%     consider using ParametricBoostrap.m for that.
%   tail: [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference)
%   glmeArgs: cell array of arguments to pass into fitglme (as glmeArgs{:})
%       e.g. {'link','logit','distribution','binomial','DummyVarCoding','reference'}
%   chan_hood = [nCh * nCh] symmetric binary matrix indicating which channels are
%               neighbors. If chan_hood(a,b)=1, then Channel A and Channel
%               B are nieghbors: spatial_neighbors(chanLocs, 0.61, [])
%   nPerms = number of permutations to run, default=5000
%   useCluster: 1=cluster-correct (Default), 0= instead just do normal
%       permutation-testing correction
%
% Outputs:
%   corrP: [coefficients, time-samples, chans] matrix of clustered permutation p-values
%   trueT: t-values for the 'true' effect, i.e. not permuted
%   trueB: beta coefficients for the same
%   trueP: p-values from 'true' regressions (i.e. not permuted or corrected)
%   permP: [coeff, times, chans, perms] p-values from each permutation test
%   permT: t-values from each permutation
%   permB: beta-coefficients from each permutation
%
% John Grogan, 2023

%% check inputs

if ~ismatrix(dvMat) || ndims(dvMat)>3
    error('dvMat must be 2D or 3D matrix');
end
if ~exist('useCluster','var') || isempty(useCluster)
    useCluster = 0;
end
if ndims(dvMat)==3 && useCluster &&  (~exist('chanLocs','var')  || isempty(chan_hood))
    error('dvMat has 3rd dimension for channel, but chanlocs is not given');
end

[n, nT, nCh] = size(dvMat);

if height(dataTab)~=n
    error('dvMat and dataTab have different numbers of rows');
end

if useCluster && (~exist('chan_hood','var') || isempty(chan_hood))
    if nCh==1
        chan_hood = false;
    else
        error('chan_hood must be given: see FindClustersLikeGND.m');
    end
end


%% remove NaN-rows to speed things up
% check they match first

nanRows = all(isnan(dvMat),[2 3]);
nanRowsTab = all(isnan(table2array(dataTab)),2);
if ~all(nanRows==nanRowsTab)
    warning('nan rows in dvMat and dataTab do not match up');
    keyboard;
end

dvMat(nanRows,:,:) = [];
dataTab(nanRows,:) = [];
% n = size(dvMat,1); % update

% also remove unused columns from table?

%% run one to check how many coefficients are returned

dataTab.DV = dvMat(:,1,1);
f = fitglme(dataTab, formula, glmeArgs{:});
coeffNames = f.CoefficientNames;

coeffInds = ismember(coeffNames, coeffToPermute);
disp(coeffNames(coeffInds));

nC = sum(coeffInds);

%% run the true regressions

[trueT, trueB, trueP] = deal(NaN(nC, nT, nCh)); % [nC, times chans]
tic;
parfor iT = 1:nT
    if mod(iT, round(nT/100))==0; disp(iT);end
    dataTab1 = dataTab;
    dv1 = dvMat(:,iT,:);
    glmeArgs1 = glmeArgs;
    for iCh = 1:nCh
        dataTab1.DV = dv1(:,1,iCh);
        f = fitglme(dataTab1, formula, glmeArgs1{:});

        trueP(:,iT,iCh) = f.Coefficients.pValue(coeffInds);
        trueT(:,iT,iCh) = f.Coefficients.tStat(coeffInds);
        trueB(:,iT,iCh) = f.Coefficients.Estimate(coeffInds);
    end
end
toc
%% now run the permutation regressions

% need to get the indices to permute along
% i.e. want to permute trials within each person only, but across any IVS
% (e.g. if we just permute within each pulse type, then any differences
% between the pulse types are preserved)

% so just get the ppID, and permute within that
colNames = {'pp1'};
tabInds = table2array(dataTab(:, ismember(dataTab.Properties.VariableNames, colNames)));

% I don't get exactly how the Shuffle() function works, so maybe just loop
% through all combinations?

[uniqs, ~, iC] = unique(tabInds, 'rows');
% iC is the index of unique rows in the table

% e.g. so for each permutation, we will shuffle all the rows where iC==1,
% then where iC==2, etc

% or can do groupMeans, and then Shuffle along 2nd dim, then put back into
% vector? It will mix in NaNs though - may be quicker to just forloop than
% try to remove them

[permT, permB, permP] = deal(NaN(nC, nT, nCh, nPerms));

tic;
parfor iP = 1:nPerms
    if mod(iP,round(nPerms/100))==0; disp(iP); end

    dvMat1 = dvMat; % copy for parfor
    dataTab1 = dataTab;
    glmeArgs1 = glmeArgs;

    % permute within pps
    iC2 = NaN(size(iC));
    for i = 1:size(uniqs,1)
        y = find(iC==i);
        iC2(iC==i) = y(randperm(length(y)));
    end

    dvMapermT = dvMat1(iC2,:,:); % permute the matrix

    for iT = 1:nT % now test each T + Ch
        for iCh = 1:nCh
            dataTab1.DV = dvMapermT(:,iT,iCh);
            f = fitglme(dataTab1, formula, glmeArgs1{:});

            permP(:,iT,iCh,iP) = f.Coefficients.pValue(coeffInds);
            permT(:,iT,iCh,iP) = f.Coefficients.tStat(coeffInds);
            permB(:,iT,iCh,iP) = f.Coefficients.Estimate(coeffInds);
        end
    end
end
toc


%% cluster-correct the t-values

if useCluster
    df = f.DFE; % degrees of freedom

    corrP = NaN(nC, nT, nCh);
    for i = 1:nC 
        % make inputs [nT nCh (nPerms)]
        corrP(i,:,:) = FindClustersLikeGND(shiftdim(trueT(i,:,:),1), shiftdim(permT(i,:,:,:),1), chan_hood, tail, df)'; %[times chans]
    end
    % do I need to adjust the t-statistic also, like permutationOLS does?
    % can check what the DMGroppe functions do

else
    % if not using clustering, just count how often it passes the max/min

    % this seems to just use the overall distribution per channel, but it
    % should probably be across time and

    %     upperP = 1 - nanmean(abs(trueT) >= permT);
    %     lowerP = 1 - nanmean(-abs(trueT) < permT);
    %
    %     corrP = (upperP + lowerP) / 2;

    % this is from permutationOLS.m in matlib
    tMax = permute(max(permT,[], [2 3]),[1,4,2,3]); % [nC nPerm]

    if tail == 0 % two-tailed
        upper_t_threshold = prctile(tMax',97.5)'; % work out cutoff values of t corresponding  [nC 1]
        lower_t_threshold = prctile(tMax', 2.5)'; % to the required family-wise error rate [nC 1]
        t_test_result = bsxfun(@gt, trueT, upper_t_threshold) ...
            - bsxfun(@lt, trueT, lower_t_threshold);       % put +1 or -1 if contrast is above or below zero [nC, nT, nCh]
        t_threshold  = cat(1,upper_t_threshold, lower_t_threshold); % T_THRESHOLD ( UPPER/LOWER, REGRESSOR )
        % For how many permutations is the maximum (across samples) permuted t-statistic
        % greater than the true t-statistic (for each sample)? If it is rarely
        % above the real t, then p is small, because t is bigger than expected by
        % chance.
        p_vals_upper = mean( bsxfun(@gt,  permute(tMax,[1 3 4 2]), trueT) , 4  ); % mean across permutations, [nC, nT, nCh]
        % How often is the permuted minimum t-statistic across samples lower than
        % the true t-statistic of each sample? If it is rarely below the real t,
        % then p is small, because t is lower than expected by chance.
        p_vals_lower = mean( bsxfun(@lt, permute(tMax,[1 3 4 2]), trueT) , 4  ); % for each contrast and sample
        p_vals       = min( p_vals_upper, p_vals_lower );  % give whichever p-value is lower. Correct for 2-tails.
        % the t_test_result will tell us which direction it it significant in.
        % The returned p-value is close to 0 if either tail is significant, but
        % must be tested with a threshold of ALPHA/2.

%         % I think we can do
%         p_vals = (p_vals_upper + p_vals_lower)/2;
% 
    else % one-tailed

        % need to flip for -1/1?


        t_threshold   = prctile(tMax',95)';         % t cutoff values for given FWER [nC 1]
        t_test_result = bsxfun(@gt, trueT, t_threshold);        % is the statistic bigger than the threshold? [nC nT nCh]
        % is the t-statistic bigger than alpha-percent of the permuted max-t
        % statistics?
        p_vals        = mean( bsxfun(@gt, permute(tMax,[1 3 4 2]), trueT) , 4  ); %[nC, nT, nCh]
        % what proportion of the permuted maximum-t statistics are above the true
        % t-statistic for each sample?
    end

    corrP = p_vals;
    

end


end
