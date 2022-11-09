function [newData, newT] = myEeglabResample(data, oldRate, newRate, t, fc, df)
% function myEeglabResample(data, oldRate, newRate, fc, df)
% Copied from pop_resample in eeglab, just the algorithm bit, without all the
% EEGlab structure stuff. It also applies resample across all trials/pps in one go
% Therefore the padding includes any NaN that are around the data - need to double
% check that this isn't giving artefacts. If so, need to rewrite this to only
% run on non-nan values and to do the padding separately
%
% Inputs:
%   data = matrix of EEG data, with time in second dimension (works up to 5-dim matrix currently)
%   oldRate = current recording frequency (Hz)
%   newRate = frequency to resample to (Hz)
% Optional inputs:
%   t = timepoints, will resample those too
%   fc         - anti-aliasing filter cutoff (pi rad / sample)
%                {default 0.9}
%   df         - anti-aliasing filter transition band width (pi rad /
%                sample) {default 0.2}
%
% Outputs:
%   newData = matrix with diff number of columns corresponding to new sampling rate
%   newT = downsampled time vector too
%   

%% inputs

if ~exist('fc','var') || isempty(fc)
    fc = 0.9;% Default cutoff frequency (pi rad / smp)
end
if ~exist('df','var') || isempty(df)
    df = 0.2;% Default transition band width (pi rad / smp)
end

if fc < 0 || fc > 1 % fc in range?
    error('Anti-aliasing filter cutoff freqeuncy out of range.')
end

% finding the best ratio
[p,q] = rat(newRate/oldRate, 1e-12); % used! AW
newRate = oldRate*p/q;


nT = size(data,2);

nT2 = ceil(nT * p/q); % new size


if nargout==2
    if ~exist('t','var') || isempty(t) 
        error('new time-vector requested as output , but not inputted');
    else
        % re-calc times
        newLims = [t(1), t(1) + (nT2-1)/newRate];
        newT = linspace(newLims(1), newLims(2), nT2);
    end
end

%% padding to avoid artifacts at the beginning and at the end
% this is unaffected by size of data on each trial

% Andreas Widmann May 5, 2011

%The resample command introduces substantial artifacts at beginning and end
%of data when raw data show DC offset (e.g. as in DC recorded continuous files)
%when MATLAB Signal Processing Toolbox is present (and MATLAB resample.m command
%is used).
%Even if this artifact is short, it is a filtered DC offset and will be carried
%into data, e.g. by later highpass filtering to a substantial amount (easily up
%to several seconds).
%The problem can be solved by padding the data at beginning and end by a DC
%constant before resampling.

%         N = 10; % Resample default
%         nPad = ceil((max(p, q) * N) / q) * q; % # datapoints to pad, round to integer multiple of q for unpadding
%         tmpeeglab = resample([data(ones(1, nPad), :); data; data(end * ones(1, nPad), :)], pnts, new_pnts);

% Conservative custom anti-aliasing FIR filter, see bug 1757
nyq = 1 / max([p q]);
fc = fc * nyq; % Anti-aliasing filter cutoff frequency
df = df * nyq; % Anti-aliasing filter transition band width
m = pop_firwsord('kaiser', 2, df, 0.002); % Anti-aliasing filter kernel
b = firws(m, fc, windows('kaiser', m + 1, 5)); % Anti-aliasing filter kernel
b = p * b; % Normalize filter kernel to inserted zeros
%         figure; freqz(b, 1, 2^14, q * 1000) % Debugging only! Sampling rate hardcoded as it is unknown in this context. Manually adjust for debugging!

% Padding, see bug 1017
nPad = ceil((m / 2) / q) * q; % Datapoints to pad, round to integer multiple of q for unpadding

%% pad now, and then run across entire matrix at once
startPad = repmat(data(:, 1, :,:,:), [1 nPad 1 1 1]);
endPad = repmat(data(:, end, :,:,:), [1 nPad 1 1 1]);

newData = [startPad, data, endPad];

newData = resample(newData, p,q,b,'dimension',2); % work across time

% Remove padding
nPad = nPad  * p / q; % # datapoints to unpad
newData = newData(:, (nPad + 1):(end - nPad), :,:,:); % Remove padded data



end
