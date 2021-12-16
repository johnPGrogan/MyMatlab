function cppSlopes = FitCPPSlope(cppTraces, slopeWindow, epochTimes)
% function FitCPPSlope(cppTraces, slopeWindow, epochTimes)
% Fit a line to the CPP within a window, return just the slope of the line
% Use regress() as it allows for NaN within trace
% 
% Inputs:
%   cppTraces: [nPP, nTimepoints, nTrials] matrix of CPP traces
%   slopeWindow: [min max] times (ms) of window to use
%   epochTimes: 1:nTimespoints vector of times (ms) for each sample
% 
% Outputs:
%   cppSlopes: [nPP, nTrials] matrix of slopes of CPP within the window
% 
% John Grogan, 2021.


% slopeWindow = [-500 -200];
slopeInds = isBetween(epochTimes,slopeWindow);
x = epochTimes(slopeInds);

[nPP, ~, nTr, nConds] = size(cppTraces);
cppSlopes = NaN(nPP,nTr,nConds,2);
for iPP = 1:nPP
    for iTr = 1:nTr
        for iCond = 1:nConds
%             cppSlopes(iPP,iTr,iCond,:) = polyfit(x, sq(cppTraces(iPP,slopeInds,iTr,iCond))', 1);
            if sum(~isnan(cppTraces(iPP,slopeInds,iTr,iCond))) > 1
                cppSlopes(iPP,iTr,iCond,:) = regress(sq(cppTraces(iPP,slopeInds,iTr,iCond)), [ones(length(x),1) x']);
            end
        end
    end
end

cppSlopes(:,:,:,1) = []; % remove intercept param


