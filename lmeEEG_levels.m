function corrP = lmeEEG_levels(eegMatrix, levelsToTest, behTab, formula, nPerms, tail, chan_hood, times, colours, yVals)
% function corrP = lmeEEG_levels(dataByFac, levelsToTest, behTab, formula, nPerms, times, colours, yVals)
% 
% Wrapper function for lmeeeg_onechan, calls that once per levelsToTest,
% allowing a range of posthoc/subsets to be run, and plots pbar() for
% significant ones.
% levelsToTest can be a cell array [1 nTests] of 3rd dimension indices to use for each
% comparison, like myPermOLS(), in which case eegMatrix should be [nPP nT nFacLevels nTrials]
% size, and behTab should be empty (nTrials can be 1 for aggregate data), OR...
% It can be a matrix [nPP*nFacLevels*nTr nTests] of row=indices to include in each one, in
% which case eegMatrix should be [nPP*nFacLevels*nTr nT nCh] matrix, and
% behTab should be table with same number of rows.
% 
% Inputs:
%   eegMatrix = EITHER:
%     1. matrix of size [nPP, nTimes, nLevels, nTr] (1 channel only), OR
%     2. matrix of size [nPP*nTr, nTimes, nChannels] - in which case, you 
%        must pass in behTab too!! (if multiple channels, will not plot
%        pbar)
%   levelsToTest = EITHER:
%     1. cell array of levels of factor to test (indices of 3rd dim of 
%        eegMatrix (if not given, will test all nLevels). [1 nTests]
%     2. matrix of row indices [nPP*nTr, nTests], one test run per column 
%        subset.
%   behTab = [OPTION 1 only] table with behavioural data and 'pp1' in, if passed
%       in, assumes that eegMatrix is a matrix or table of same height
%   formula = [OPTIONAL] regression formula with Random Effects (DV must be 'EEG'),
%     and fixed-effects either just '1 + fac' if not using behTab, or
%     variables in behTab if that is given. default is 'EEG ~ 1 + fac + (1 | pp1)'
%   nPerms = number of permutations to run per test (Same for all)
%   tail = tail of test to run, default=0 = two-tailed, -1 or 1 for
%     one-sided directional tests
%   chan_hood = if 1 channel given, is set to false, otherwise should be:
%     chan_hood = spatial_neighbors(eeg.chanlocs, 0.61, []);
% These params are given only if you want to call pbar to plot sig
% clusters, and only 1 channel given:
%   times = [optional] xValue times to use for pbar, if not given, pbar not
%     called
%   colours = [optional] colours to use for each test [nTests, 3]. default
%     is [0 0 0] and then figure colour order
%   yVals = [optional] yValues to plot pbars at [1 nTests]
% 
% Outputs:
%   corrP = [nTests, nT, nCh] pvalue outputs from lmeeeeg_onechan
% 

if exist('behTab','var') && ~isempty(behTab)
    % if this is passed in, eegMatrix should be [nPP*nL*nTr, nT, nCh] matrix
    % and behTab should be [nPP*nTr, nBehVars] table
    [nRows, nT, nCh] = size(eegMatrix);
    assert(nRows == height(behTab), 'eegMatrix and behTab have diffferent number of rows');

    assert(ismatrix(levelsToTest) && ~isempty(levelsToTest), 'levelsToTest must be matrix if behTab given');

    nTests = size(levelsToTest,2);

else % eegMatrix is 3/4D matrix, behTab is empty
    [~, nT, nL, ~] = size(eegMatrix);
    nCh = 1; % force to 1

    if ~exist('levelsToTest','var') || isempty(levelsToTest)
        levelsToTest = {1:nL}; % default to testing all levels, will use diff wave if nL<3
    end
    assert(iscell(levelsToTest), 'levelsToTest ,ust be cell array of indices if behTab not given')

    nTests = length(levelsToTest);
end

% check sizes match for pbar
if exist('times','var') && ~isempty(times)
    assert(nT == length(times), 'times and data have different lengths');
end


% use lmeeeg_onechan defaults
if ~exist('tail','var') || isempty(tail); tail = []; end
if ~exist('chan_hood','var') || isempty(chan_hood); chan_hood = []; end 

%% run each test
% assumes only 1 coefficient (excl intercept) per test

corrP = NaN(nTests, nT, nCh); % corrected pvluaes
for iT = 1:nTests
    if iscell(levelsToTest)
        % use matrix version, no behTab passed in
        corrP(iT,:) = lmeEEG_oneChan(eegMatrix(:,:,levelsToTest{iT},:), [],...
            formula, nPerms, tail, chan_hood);
    
    else % use rowIndices version
        corrP(iT,:) = lmeEEG_oneChan(eegMatrix(levelsToTest(:,iT),:),...
            behTab(levelsToTest(:,iT),:), formula, nPerms, tail, chan_hood);

    end
end

%% 
if exist('times','var') && ~isempty(times)
     % plot pbars

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
end