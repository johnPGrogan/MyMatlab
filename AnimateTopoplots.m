function AnimateTopoplots(data, times, timeInds, chanlocs, varNames, mapLims, saveName, cmap, cb, imwriteArgs)
% function AnimateTopoplots(data, times, timeInds, chanlocs, varNames, mapLims, saveName, cmap, imwrteArgs)
% Draw multiple animated topoplots on one figure (e.g. one subplot per variable)
% 
% Inputs:
%   data = [chans times OTHER] data to plot. will do one subplot per OTHER
%   times = vector of time points
%   timeInds = vector of indices of times to plot (one frame per index)
%   chanLocs = chanlocs info from ERP/EEG
% Optional:
%   varNames = names of variables, to match number of OTHERs
%   mapLims = [min max] colour map limits. default is symmetrical minmax
%   saveName = 'name.gif' including path, to save gif as. './test.gif' is
%       default
%   cmap = colourmap e.g. 'jet' or crameri('vik'). default = jet
%   cb = 1 = draw colour bar on final subplot. default = 1
%   imwriteArgs = cell containing param-val pairs to pass into imwrite().
%       default is {'DelayTime',0.05,'LoopCount',inf}

if ~exist('varNames','var') || isempty(varNames)
    varNames = repmat({''}, size(data,3), 1);
end
if ~exist('mapLims','var') || isempty(mapLims)
    mapLims = [-1 1] .* max(abs(data(:,timeInds,:,:,:,:)),[],'all');
end
if ~exist('saveName','var') || isempty(saveName)
    saveName = './test.gif';
end
if ~exist('cmap', 'var') || isempty(cmap)
    cmap = 'jet';
end
if ~exist('cb','var') || isempty(cb)
    cb = 1;
end
if ~exist('imwriteArgs','var') || isempty(imwriteArgs)
    imwriteArgs = {'DelayTime',0.05,'LoopCount',inf};
end
%%
    
hsig = figure();
set(hsig,'Position', get(0,'Screensize'));
% ERP.bindata = permute(nanmean(varCorr(:,:,1,1,:,:),1),[6,2,5,1,3,4]); %[chans times vars], average over pps. 1 cond + drug
ERP.bindata = data;
ERP.chanlocs = chanlocs;
ERP.times = times;
nVars = size(data,3);


% timeInds = 1:1:666;
fra = 1;
for iT = 1:length(timeInds)
    clf
    SubTopoPlot(ERP, timeInds, iT, varNames, cmap, nVars, mapLims, cb);
    
    % save
    if fra==1
        f = getframe(hsig);
        pause(0.01)
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        pause(0.01)
        im(1,1,1,iT) = 0;
    end
    
    fra = fra + 1;
    f   = getframe(hsig);
    pause(0.01)
    im(:,:,1,fra) = rgb2ind(f.cdata, map,'nodither');
    pause(0.02)
    
end

pause(0.1)
%
imwrite(im, map, saveName, imwriteArgs{:});% 


end

function SubTopoPlot(ERP, timeInds, iT, varNames, cmap, nVars, mapLims, cb)
% Draw multiple topoplots on one figure
    [r,c] = GetSubPlotShape(nVars);
    for i = 1:nVars
        subplot(r,c,i)
        topoplot( ERP.bindata(:,timeInds(iT),i), ERP.chanlocs,...
            'style', 'both', 'plotrad',.55, 'headrad', .5,'emarker', {'.','k',[],1},...
            'numcontour', 6, 'maplimits', mapLims, 'colormap', cmap,'electrodes', 'on', 'nosedir', '+X');
        
        title(varNames{i});
    end
    if cb
        colorbar;
    end
    SuperTitle(sprintf('t = %.0f ms', ERP.times(timeInds(iT))),{'Position',[.33 1 .5]});
    
end