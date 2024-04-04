function corrP = lmeEEG_categorical_pbar(eegMatrix, behTab, formula, nPerms, tail, chan_hood, categOrder, contrasts, times, colours, yVals, pbarsToPlot)
% function p = lmeEEG_categorical_pbar(eegMatrix, behTab, formula, nPerms, tail, categOrder, contrasts, times, colours, yVals, pbarsToPlot)
% 
% Wrapper function for lmeEEG_cateogircal that also does pbar plotting
% As categorical predictor is used, will give pbar per level of category
% (except reference), plus any contrasts tested also.
% 
% First 8 inputs are the same as lmeEEG_categorical, so see help for that
% and the final 3 inputs are below (if not given, no pbar is plotted, so is
% the same as calling lmeEEG_categorical.
% 
%     times: [1 nTimes] vector of times to plot pbar at. if not given, no
%       pbars will be plotted
%     colours = colours to plot each pbar, either [n 3] rbg matrix, or [n 1]
%       vector of colour codes e.g. ['r';'b']. default is 
%        [0 0 0; get(gca, 'ColorOrder')].
%     yVals = [n 1] vector of yValues to plot each pbar at - only significant
%       pbars are plotted and non-sig are skipped, so yvals are used in order
%       of significant pbars, not order of tests given. default is min(ylim)+i
%     pbarsToPlot = indices of which pbars to plot, e.g. can just plot
%       overall contrast, or change order (to match those in colours &
%       yVals)
% 
% Outputs:
%     p = cluster-corrected pvalues [nTests, nTimes]
% 

if ~exist('times','var') || isempty(times)
    warning('no times are given, so no pbars will be plotted')
else
    assert(size(eegMatrix,2) == length(times), 'times and data have different lengths');
end


% set default inputs as empty, will use defaults in lmeEEG_categorical
if ~exist('behTab','var') || isempty(behTab); behTab = []; end
if ~exist('formula','var') || isempty(formula); formula= []; end
if ~exist('nPerms','var') || isempty(nPerms); nPerms = []; end
if ~exist('tail','var') || isempty(tail); tail = []; end
if ~exist('chan_hood','var') || isempty(chan_hood); chan_hood = []; end 
if ~exist('categOrder','var') || isempty(categOrder); categOrder= []; end
if ~exist('contrasts','var') || isempty(contrasts); contrasts = {}; end 

%% call lmeEEG_categorical


corrP = lmeEEG_categorical(eegMatrix, behTab, formula, nPerms, tail, chan_hood, categOrder, contrasts);


%% plot pbars

if exist('times','var') && ~isempty(times)
     % plot pbars
     % only some?

     if exist('pbarsToPlot','var') && ~isempty(pbarsToPlot)
         corrPOrig = corrP; % store
         corrP = corrP(pbarsToPlot,:); % discard or rearrange
     end

    isSig = any(corrP <= .05, 2);  % only plot sig ones
    if any(isSig) % only plot if there are any
        if ~exist('yVals','var') || isempty(yVals)
            yVals = min(ylim) + (1:sum(isSig)) - 1; % default yVals
        end
    
        if ~exist('colours','var') || isempty(colours)
            colours = [0 0 0; get(gca,'ColorOrder')];
        end
    
    
        f = find(isSig);% can skip out any missing ones?
        for iT = 1:length(f)
            hold on;
            pbar(corrP(f(iT),:), 'xVals',times, 'yVal', yVals(iT), 'plotargs',{'LineWidth',5,'Color',colours(f(iT),:)});
        end 
    end

    if exist('corrPOrig','var') % reset to orig order
        corrP = corrPOrig;
    end
end


