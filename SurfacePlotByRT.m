function SurfacePlotByRT(erp, rt, erpTimes, tLims, smoothWin, cLims)
% function SurfacePlotByRT(erp, rt, erpTimes, tLims, smoothWin)
% 
% Surface plot of erp sorted by RT, with smoothing applied
% Inputs:
%   erp = [nPP nT nTr] ERP, will be turned into [nPP*nTr, nT] and sorted
%    and smoothed by RT across all pps/trials
%   rt = [nPP nTr] RT, will sorted across all pps/trials
%   erpTimes = [1 nT] time-points of ERP 
%   tLims = [min max] ERP times to plot
%   smoothWin = number of trials to smooth sorted RT + ERPs by
%     (default=100)
%   cLims = imagesc colour limits
% 

%% setup

if ~exist('smoothWin', 'var') || isempty(smoothWin)
    smoothWin = 100;
end
if ~exist('cLims', 'var') || isempty(cLims)
    cLims = [];
end

[nPP, nT, nTr] = size(erp);
[~, nT2] = size(erpTimes);
[nPP2, nTr2] = size(rt);

% check sizes match
assert( nPP==nPP2 &&  nTr==nTr2 , 'size of erp + rt do not match')
assert(nT==nT2, 'erp and erpTimes have different lengths');

% sort x-limits + x-ticks
% imagesc plots by indices, so need to use xticklabels to set them

xInds = isBetween(erpTimes, tLims); % get indices to plot
x = erpTimes(xInds); % get times plotting
[~, xTickInds] = min( abs(x - (tLims(1):200:max(x))') ,[],2); % find every 200ms indices
xLabels = round(x(xTickInds),2,'significant'); % round them

% find offset of 1st point from zero, to offset RT by
x0 = -x(1)./unique(diff(erpTimes));
x0 = x0(1); % in cast there are slight rounding errors

%% reshape both - collapse across pp + tr

rt = reshape(rt, nPP*nTr, 1); %[nPP*nTr, 1]
erp = reshape(permute(erp,[1,3,2]), nPP*nTr, nT); %[nPP*nTr, nT]

%% remove NaNs

toRemove = isnan(rt) | all(isnan(erp),2);

rt(toRemove) = [];
erp(toRemove,:) = [];

%% sort by RT

[rt, latOrder] = sort(rt);
erp = erp(latOrder,:);

%% smooth

rt = movmean(rt, smoothWin, 1);
erp = movmean(erp, smoothWin, 1);

%% plot

imagesc(erp(:,xInds), cLims); % plot erp
% plot RT line
hold on;
plot(rt + x0, 1:numel(rt), '-k','LineWidth',2);
       

set(gca,'YDir','normal'); % ascending y axis 
colorbar; % show this

% change xticklables
xticks(xTickInds);
xticklabels(xLabels);

xline(x0,':'); % show t=0

ylabel('trials');


