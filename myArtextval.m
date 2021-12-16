function isArtefact = myArtextval(erp, thresh)
% copied from pop_artexval.m (eeglab) - just the main algorithm
% Marking epochs containing activity above an upper threshold and below a lower threshold
% Inputs
%   erp = [chans times trials] matrix of voltages - will be applied across
%          all channels and times, so trim before inputting
%    thresh = [min max] threshold for detections (uses <= and >=)
%
% Outputs:
%     isArtefact = logical [chans 1 trials]
%

[nCh, ~, nTr] = size(erp);


isArtefact = false(nCh,1,nTr);

for iCh=1:nCh
    
%     fprintf('%g ',iCh);
    
    for iTr = 1:nTr
        dataline = erp(iCh, : ,iTr);
        criteria1 = max(dataline)>= thresh(2);
        criteria2 = min(dataline)<= thresh(1);
        if criteria1 || criteria2
            isArtefact(iCh,1,iTr) = true;
        end
    end
end