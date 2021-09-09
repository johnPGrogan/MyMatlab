function w = wITPC(data, dims, weights)
% function w = wITPC(data, dims, weights);
% Calculate weighted ITPC, across any dims. If weights are not given, is
% the same as ITPC
% Inputs:
%   data = matrix of complex time-freq signal
%   dims = dimensions to take average over, can be scalar or vector
%   weights = broadcastable matrix/vector of weights to multiply the
%       lengths of the vectors by
% Outputs:
%   w = squeezed matrix (after the averaging) of wITPC. abs(w) will give
%       coherence, and angle(w) will give coherent angle
% 
% John Grogan, 2021.

if ~exist('dims','var') || isempty(dims) || any(dims > ndims(data))
    error('must supply dims that match ndims(data)');
end
if ~exist('weights','var') || isempty(weights)
    weights = 1;
end

% eik = data ./ sqrt(data.* conj(data)); % get polar representation of phase angles
% breik = eik .* weights; % multiply by weights (only changes vector magnitude)
% w = squeeze(nanmean(breik, dims)); % average 

% witpc = sq( nansum(hilw ./ sqrt(hilw .* conj(hilw)), dim) / size(x,dim) );

% w1 = sq(nanmean( exp(1i * angle(data)) .* weights, dims));


eik2 = exp(1i * angle(data)); % polar representation of phase angles
breik2 = eik2 .* weights; % multiply by weights (changes magnitude not angle)
w = sq(nanmean(breik2, dims)); % average them

% equals(round(w,4), round(w1,4))

% max(abs(w) - abs(w1), [], 'all')
% max(angle(w) - angle(w1), [], 'all')
end