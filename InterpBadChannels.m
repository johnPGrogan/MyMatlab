function InterpBadChannels(dataFile, interpFolder)

doPlots = 0;
outFolder = './Saves';

if ~exist('dataFile','var') || isempty(dataFile)
    dataFile = fullfile(outFolder, 'ExtractEpochs.mat');
end
load(dataFile);

if ~exist('interpFolder','var') || isempty(interpFolder)
    interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Interp';
end
% if ~isfield(fileInfo, 'interpFolder')
    fileInfo.interpFolder = interpFolder;
    save(dataFile, '-append','fileInfo');
% end

% These came from running ArtifactCheck on the *raw.mats:
load('ChansToInterp2.mat'); %
for iPP = 1:fileInfo.nPP

if ~exist(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} '_intp.mat']), 'file') % don't re-do
    disp(iPP)
    
    doInterp(iPP, chansToInterp, fileInfo.rawFolder, fileInfo.ppID, fileInfo.interpFolder, eeg.nChans, eeg.nSamples, eeg.chanlocs, isBetween(eeg.epochTimes, eeg.blWin));
else
    fprintf('\n %d already interpolated', iPP);
end

end

disp('done!');

end

function doInterp(iPP, chansToInterp, rawFolder, ppID, interpFolder, nChans, nSamples, chanlocs, blInds)

    disp('loading...');
    epoch = load(fullfile(rawFolder, [ppID{iPP} '_raw']));     

    
    % apply across entire dataset
    badChans = chansToInterp.(ppID{iPP});   % get the bad channels recorded during ArtifactCheck for that block
    if ~isempty(badChans) % only do this if bad chans, otherwise just rereference
        

        load('blankEEG.mat');  % load a clean-slate EEG structure that has no data, but has all the right fields for EEGLAB to jazz with it
        EEG.nbchan = nChans; % set number of electrodes to just the number of scalp electrodes - we will not interpolate based on the externals!
        EEG.data = epoch.erp(1:nChans,:,:); % feed it the epoched data for this block, again only the scalp channels in the cap
        EEG.pnts = nSamples; % it seems to need this too - epoch length in sample pts
        EEG.trials = size(EEG.data,3); % number of trials
        EEG.chanlocs = chanlocs(1:nChans); % it needs channel locations too

        % un-baseline for now (will mean that baseline period gets interpolated too now)
        EEG.data = EEG.data + epoch.baseline(1:nChans,:,:);

        disp('interpolating...');
        EEG=eeg_interp(EEG,badChans,'spherical'); % this line does the actual interpolation

        % do I also need to interpolate the baselines - used in CSD

        if 1% doPlots == 1
            clf;subplot(1,2,1);errorBarPlot(permute(epoch.erp(badChans,:,:),[3 2 1]), 'area',1); title('raw'); 
            ylim([-30 30]);
            % this was plotting un-baselined, so often looked worse
%             subplot(1,2,2);errorBarPlot(permute(EEG.data(badChans,:,:),[3 2 1]), 'area',1);title('interpolated');
%     %             makeSubplotScalesEqual(1,2);
%             ylim([-30 30]);
%             SuperTitle(iPP);
%             drawnow;
        end
        epoch.erp(1:nChans,:,:) = EEG.data; % now replace the relavant parts of the big 'erp' matrix with the interpolated version
            % Note the externals will still be sitting there in channels 129-136, unaltered.
    
        % re-baseline
        epoch.baseline(1:nChans,:,:) = nanmean(EEG.data(:, blInds, :),2);
        epoch.erp(1:nChans,:,:) = EEG.data - epoch.baseline(1:nChans,:,:); % subtract baseline
        
        % show re-baselined data 
        if 1 % doPlots==1
            subplot(1,2,2);errorBarPlot(permute(epoch.erp(badChans,:,:),[3 2 1]), 'area',1); title('raw'); 
    %             makeSubplotScalesEqual(1,2);
            ylim([-30 30]);
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

