function pval = FindClustersLikeGND_F(trueFVals, permFVals, chan_hood, tail, df)
% Function to find significant clusters after having already run a set of
% permutation tests, using the resulting null F-statistic distribution.
% Since F-statistics are >0, needed to adjust the t-statistic stuff to use
% different thresholds etc.
% The code below is copied from the Mass Univariate Toolbox (DMGroppe)
% functions clustGND.m & clust_perm1, and requires that toolbox be on the path.
%
% Inputs:
%   trueFVals = the F-values for the 'true' effect, i.e. without any
%       permutation, in [nTimes, nChans] size
%   permFVals = the F-values from the permutation tests [nTimes, nChans, nPerms]
%   chan_hood = [nCh * nCh] symmetric binary matrix indicating which channels are
%               neighbors. If chan_hood(a,b)=1, then Channel A and Channel
%               B are nieghbors: spatial_neighbors(chanLocs, 0.61, [])
%   tail = ignored, legacy argument to match FindClustersLikeGND.m. Since
%     F-distribution is >0, forces a one-tailed test
%   df = degrees of freedom to use [df1, df2]
%
%
% Outputs:
%   pval = [chans times] matrix of clustered permutation pvals
% 
% John Grogan, 2022.

tail = 1; % force to this

assert(~(any(trueFVals < 0, 'all') | any(permFVals<0, 'all')), 'F values cannot be negative');


%% get info

[nTimes, nChans, nPerms] = size(permFVals);

%% get channel neighbours

% chan_dist = 0.61; % default from the toolboxes,
% head_radius = [];
% 
% chan_hood = spatial_neighbors(chanLocs(1:nChans), chan_dist, head_radius);

%% set desired f_thresh

% tail = 0;

thresh_p = .05; % desired alpha
fwer = .05;

% df = 29;%18576; % df in fitglme (small differences don't matter)

if tail
    %one tailed test
    thresh_f = finv(1-thresh_p, df(1), df(2)); % 1-alpha to get upper number
% else
%     %two tailed test
%     thresh_t=tinv(thresh_p/2,df);
end


%% on each perm, get cluster

% nPerms = 100;
mn_clust_mass = zeros(1,nPerms);

parfor iPerm = 1:nPerms

    f = permFVals(:,:,iPerm)'; %[chans times] % get the permuted f-vals

    % find clusters of positive values
    [clust_ids, n_clust] = find_clusters(f, thresh_f, chan_hood, 1); % find the clusters that pass thresh


    %get most extremely positive F-score 
    mn_clust_mass(iPerm) = find_mn_mass(clust_ids, f, n_clust);

end

% Estimate true FWER of test
%one tailed
fmx_ptile=prctile(mn_clust_mass,100*(1-fwer));
est_alpha=mean(mn_clust_mass >= fmx_ptile);

fprintf('Desired family-wise error rate: %f\n',fwer);
fprintf('Estimated actual family-wise error rate: %f\n',est_alpha);

%% get pvals

trueFVals = trueFVals'; % make into [chans times] for this bit

pval=ones(nChans,nTimes); % will rotate back afterwards

%upper tailed
[clust_ids, n_clust] = find_clusters(trueFVals, thresh_f, chan_hood, 1); % find clusters above this
clust_info.pos_clust_pval=ones(1,n_clust);
clust_info.pos_clust_mass=zeros(1,n_clust);
clust_info.pos_clust_ids=clust_ids;
for a=1:n_clust
    use_ids=find(clust_ids==a);
    clust_mass=sum(trueFVals(use_ids));
    clust_p=mean(mn_clust_mass >= (clust_mass)); % p(>=)
    pval(use_ids)=clust_p;
    clust_info.pos_clust_pval(a)=clust_p;
    clust_info.pos_clust_mass(a)=clust_mass;
end

pval = pval';% rotate back into [nTimes nChans]


end



%% subfunc

function mn_clust_mass=find_mn_mass(clust_ids,data_t,n_clust)
% find positive clusters
% copied from clust_perm1.m in Mass_Univariate_Toolbox

mn_clust_mass=0;

%looking for most positive cluster mass
for z=1:n_clust
    use_ids=(clust_ids==z);
    use_mass=sum(data_t(use_ids));
    if use_mass > mn_clust_mass
        mn_clust_mass=use_mass;
    end
end

end