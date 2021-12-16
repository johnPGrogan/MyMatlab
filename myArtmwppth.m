function isArtefact = myArtmwppth(erp, thresh, windowSize, stepSize, fs)
% copied from pop_artmwppth (eeglab)
% Marking epochs containing peak to peak activity, within a moving window,
% that is greater than a given threshold.
%
% Inputs:
%   erp = [chans time tr] matrix. will run on all chans and times, so trim
%       before inputting
%   thresh = threshold of difference in peaks within window
%   windowSize = time of window to use (ms)
%   stepSize = step size (ms)
%   fs = sampling rate in Hz e.g. 512
%
% Outputs:
%   isArtefact = true when artefact detected [chan 1 tr]
%


%%

[nCh, nT, nTr] = size(erp);
winSamples  = floor(windowSize*fs/1000); % get sizes in samples
stepSamples  = floor(stepSize*fs/1000);

%%

isArtefact = false(nCh,1,nTr);
for iCh = 1:nCh
%     fprintf('%g ',iCh);
    for i = 1:nTr
        for j = 1:stepSamples:nT-(winSamples-1)
            w1  = erp(iCh, j:j+winSamples-1,i);
            vs  = abs(max(w1)- min(w1));
            
            if vs > thresh 
                isArtefact(iCh,1,i) = true;
                break;
            end
        end
    end
end
