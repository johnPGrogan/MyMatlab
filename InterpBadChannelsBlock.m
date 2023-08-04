function InterpBadChannelsBlock(dataFile, interpFolder)
% channel interpolation done within each block
% Inputs:
%   dataFile = name of file with fileInfo + eeg saved in it
%   interpFolder = path of folder to save interpolated data in


if ~exist('dataFile','var') || isempty(dataFile)
    outFolder = './Saves/';
    dataFile = fullfile(  outFolder, 'ExtractEpochsDelayedConf.mat');
%     error('no fileInfo file given');
end
load(dataFile);

if ~exist('interpFolder','var') || isempty(interpFolder)
%     error('no interp folder given');
    interpFolder = '../../../../../OneDrive - TCDUD.onmicrosoft.com/ProjectsOD/DelayedConf/Analysis/Data/Interp/';
end
% if ~isfield(fileInfo, 'interpFolder')
    fileInfo.interpFolder = interpFolder;
    save(dataFile, '-append','fileInfo'); % update this
% end

% These came from running ArtifactCheck on the *raw.mats:
load('ChansToInterpBlock.mat');
for iPP = 1:fileInfo.nPP

if 1%~exist(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} '_intp.mat']), 'file') % don't re-do
    disp(iPP)
    
    doInterp(iPP, chansToInterp, fileInfo.rawFolder, fileInfo.ppID, fileInfo.interpFolder, eeg.nChans, eeg.nSamples, eeg.chanlocs, isBetween(eeg.epochTimes, eeg.blWin));
else
    fprintf('\n %d already interpolated', iPP);
end

end

disp('done!');

end

function doInterp(iPP, chansToInterp, rawFolder, ppID, interpFolder, nChans, nSamples, chanlocs, blInds)

    doPlots = 1;
    disp(['loading ' [ppID{iPP} '_raw']]);
    epoch = load(fullfile(rawFolder, [ppID{iPP} '_raw']),'erp','baseline','ch2interp','blockNum');     

    badChans = chansToInterp.(ppID{iPP});   % get the bad channels recorded during ArtifactCheck for that block

    % check ch2interp match up
    if ~equals(epoch.ch2interp, badChans)
        error('bad-channels in chansToInterp does not match those in the _raw data file');
    end
    
    if ~isempty(badChans) % only do this if bad chans, otherwise just rereference
        

        load('blankEEG.mat');  % load a clean-slate EEG structure that has no data, but has all the right fields for EEGLAB to jazz with it
        EEG.nbchan = nChans; % set number of electrodes to just the number of scalp electrodes - we will not interpolate based on the externals!
        EEG.data = epoch.erp(1:nChans,:,:); % feed it the epoched data for this block, again only the scalp channels in the cap
        EEG.pnts = nSamples; % it seems to need this too - epoch length in sample pts
        EEG.trials = size(EEG.data,3); % number of trials
        EEG.chanlocs = chanlocs(1:nChans); % it needs channel locations too


        % un-baseline for now (will mean that baseline period gets interpolated too now)
        EEG.data = EEG.data + epoch.baseline(1:nChans,:,:);

        disp('interpolating per block:');
        nBlocks = max(epoch.blockNum);
        for i = 1:nBlocks
            fprintf(' %d, ', i);
            blockInds = epoch.blockNum == i; % get trials in block
            EEGBlock = EEG; % copy
            EEGBlock.data = EEG.data(:,:,blockInds); % split
            EEGBlock.trials = size(EEGBlock.data,3); % number of trials

            % run interpolation
            EEGBlock = eeg_interp(EEGBlock,badChans{i},'spherical'); % this line does the actual interpolation

            EEG.data(:,:,blockInds) = EEGBlock.data; % put into main structure
        end

        % do I also need to interpolate the baselines - used in CSD

        if doPlots == 1
            clf;
            for i = 1:nBlocks
                if isempty(badChans{i}); continue; end % skip
                subplot(nBlocks,2,i*2-1); % left column
                blockInds = epoch.blockNum == i; % get trials in block
                h = errorBarPlot(permute(epoch.erp(badChans{i},:,blockInds),[3 2 1]), 'area',1); 
                title(sprintf('raw, block %d', i)); 
                ylim([-50 50]);
            end
        end
        epoch.erp(1:nChans,:,:) = EEG.data; % now replace the relavant parts of the big 'erp' matrix with the interpolated version
            % Note the externals will still be sitting there in channels 129-136, unaltered.
    
        % re-baseline
        epoch.baseline(1:nChans,:,:) = nanmean(EEG.data(:, blInds, :),2);
        epoch.erp(1:nChans,:,:) = EEG.data - epoch.baseline(1:nChans,:,:); % subtract baseline
        
        % show re-baselined data 
        if doPlots==1
            for i = 1:nBlocks
                if isempty(badChans{i}); continue; end % skip
                subplot(nBlocks,2,i*2); % left column
                blockInds = epoch.blockNum == i; % get trials in block
                h = errorBarPlot(permute(epoch.erp(badChans{i},:,blockInds),[3 2 1]), 'area',1); 
                title(sprintf('interpolated, block %d', i)); 
                ylim([-50 50]);
            end
            SuperTitle(iPP);
            drawnow;
        end
    end
    
    % average reference
    epoch.ref = nanmean(epoch.erp,1);
    epoch.erp = epoch.erp - epoch.ref;
    
    epoch.chansInterpolated = badChans; % store this
    
    % now re-save inerpolated version:
    disp('saving...');
    save(fullfile(interpFolder, [ppID{iPP} '_intp']),...
        '-v7.3', '-struct', 'epoch');
    
end

