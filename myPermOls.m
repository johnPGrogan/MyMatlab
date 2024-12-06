function p = myPermOls(dataByFac, levelsToTest, permArgs, times, colours, yVals, plotArgs)
% function myPermOls(dataByFac, levelsToTest, permArgs, times, colours, yVals, plotArgs)
% 
% Wrapper function for permutationOLS, allows easy running of multiple
% difference waves or comparisons, and plotting of the pbars
% 
% Inputs:
%   dataByFac = [nPP, nT, factor (nLevels)] matrix of time-series data.
%   levelsToTest = cell array of levels of factor to test, if not given,
%     will test all nLevels. If 2 levels given, will calculate difference
%     wave using diff(x, [], 3). If 1 level, assumes it is already a
%     difference wave. If >2 levels, will create a design matrix with the
%     levels and test a linear contrast of that factor.
%   permArgs = cell array of args to pass to permutationOLS(), given as
%     permArgs{:}, 
%   times = [optional] xValue times to use for pbar, if not given, pbar not
%     called
%   colours = [optional] colours to use for each test [nTests, 3]. default
%     is [0 0 0] and then figure colour order
%   yVals = [optional] yValues to plot pbars at [1 nTests]
%   plotArgs = cell array of pbar arguments (apart from color). default is
%     {'LineWidth',5}
% 
% Outputs:
%   p = [nTests, nT] pvalue outputs from permutaitonOLS
% 

[nPP, nT, nL] = size(dataByFac);


if ~exist('levelsToTest','var') || isempty(levelsToTest)
    levelsToTest = {1:nL}; % default to testing all levels, will use diff wave if nL<3
end
nTests = length(levelsToTest);

if ~exist('permArgs','var') || isempty(permArgs)
    permArgs = {};
end


p = NaN(nTests, nT);
for iT = 1:nTests
    if length(levelsToTest{iT})>2
        disp('multiple level factor: '); disp(levelsToTest{iT});
        % passing multiple level factor in
        nL = length(levelsToTest{iT}); % update

        d = reshape(permute(dataByFac(:,:,levelsToTest{iT}),[1,3,2]),nPP*nL,nT); % make [nPP*nL, nT]
        DES = [flat(ones(nPP,nL))  zscore(flat(repmat(1:nL,[nPP,1]))) ]; % design matrix [1, factor]
        group = repmat( (1:nPP)',nL,1); % pp index, will permute factor within pps
        contrasts = [0 1]; % ignore pp, just test factor
        [ ~, p(iT,:)] = permutationOLS(d, DES, contrasts, group, permArgs{:});  

    elseif length(levelsToTest{iT})==2
        disp('creating difference wave'); disp(levelsToTest{iT});
        % make difference wave
        d = diff(dataByFac(:,:,levelsToTest{iT}),[],3); %[nPP nT]

        % will randomly flip +1/-1 within pps
        [~, p(iT,:)] = permutationOLS(d, [],[],[], permArgs{:}); 


    else % assuming already a difference wave
        disp('assuming a difference wave'); disp(levelsToTest{iT});
        d = dataByFac(:,:,levelsToTest{iT});
        [~, p(iT,:)] = permutationOLS(d, [],[],[], permArgs{:}); 

%     else
%         error('levelsToTest is wrong');
    end
end

%% 
if exist('times','var') && ~isempty(times)
     % plot pbars

     if ~exist('plotArgs','var') || isempty(plotArgs)
         plotArgs = {'LineWidth',5}; % colour will be set below
     end

    isSig = any(p < .05, 2);  % only plot sig ones
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
            pbar(p(f(iT),:), 'xVals',times, 'yVal', yVals(iT), 'plotargs',{'Color',colours(f(iT),:), plotArgs{:}});
        end 
    end
end