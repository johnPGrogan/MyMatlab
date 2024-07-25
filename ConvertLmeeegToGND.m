function GND = ConvertLmeeegToGND(betas, t_obs, corrP, se, times, chanNames, chanlocs, bin_names, exp_name)
% turns outputs from lmeeeg into GND structure that can be passed into
% gui_erp. betas/t_obs,corrP,se should all be size [nCoeffs, nTimes, nChans]
% and it will create a 'bin' for each coefficient. It uses beta-coeffs in
% place of voltage, and shows t-scores on bototm panel. 
% You can then plot this via DMGroppe mass univariate toolbox:
% gui_erp('initialize','GNDorGRP',GND,'t_test',coeffIndex,'stat','t','verblevel',1);
% 
% Inputs:
%   betas = beta-coefficients output from lmeeeg_OneChan [nCoeffs, nTimes, nChans]
%   t_obs = observed t_stats from lmeeeg, same shape as beta
%   corrP = permutation-corrected p-values, same shape as beta
%   se = beta-coefficient standard errors, same shape as beat
%   times = vector of times [1 nTimes]
%   chanNames = cell array of channel names
%   chanlocs = structure array of channel locations
%   bin_names = names of each coefficient, which will be used as bin-labels
%   exp_name = [optional] string to put as exp name, used to figure title
%       (default is empty string)
% 
% Outputs:
%   GND: structure for use in gui_erp (DMGroppe mass univariate toolbox)
%   

[nCoeffs, nTimes, nChans] = size(betas);
assert(all(size(betas) == size(t_obs),'all'), 't_obs and betas are different sizes');
assert(all(size(betas) == size(corrP),'all'), 'corrP and betas are different sizes');
assert(all(size(betas) == size(se),'all'), 'se and betas are different sizes');
assert(length(times) == nTimes, 'number of timepoints in betas does not match times');
assert(length(chanNames) == nChans, 'chanNames length does not match number of channels in betas');
assert(length(chanlocs) == nChans, 'chanlocs length does not match number of channels in betas');

if ~exist('exp_name','var'); exp_name = ''; end


%% put data into structure

GND = struct();
GND.exp_desc = exp_name;
GND.chanlocs = chanlocs;
GND.time_pts = times;

% put data in, as [nChans, nTimes, nCoeffs]
GND.grands = permute(betas,[3,2,1]); % use beta-coeffs as 'voltage'
GND.grands_t = permute(t_obs,[3 2 1]); % t-stats
GND.grands_stder = permute(se,[3 2 1]); 

% add in p-value stuff, per coeff
for i = 1:nCoeffs
    GND.bin_info(i).condcode = i;
    GND.bin_info(i).bindesc = sprintf('%s', bin_names{i});
    
    GND.t_tests(i).bin = i;
    GND.t_tests(i).time_wind = times([1 end]);
    GND.t_tests(i).used_tpt_ids = 1:nTimes;
    GND.t_tests(i).used_chan_ids = 1:nChans;
    GND.t_tests(i).include_chans = reshape(chanNames,[],1); % as column
    GND.t_tests(i).null_mean = 0; % placeholder
    GND.t_tests(i).mean_wind = 'no'; % placeholder
    GND.t_tests(i).crit_t = NaN; % use p-vals to calculate it
    GND.t_tests(i).adj_pval = squeeze(corrP(i,:,:))'; %[chan time]
    GND.t_tests(i).estimated_alpha = .05;
    GND.t_tests(i).desired_alphaORq = .05;
end



