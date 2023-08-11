function CSDTransform(interpFolder, csdFolder, rebaseline, dataFile)
% applies CSD transform to all files within a folder, saves them in \CSD


if ~exist('interpFolder','var') || isempty(interpFolder)
    interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Interp\';
end
disp(interpFolder);

if ~exist('dataFile','var') || isempty(dataFile)
    dataFile = '.\Saves\ExtractEpochs.mat';
end
load(dataFile,'eeg'); % get info

% load matrices that CSD function needs to run for the specific biosemi 128 cap:
load('CSD_coords'); % contains G, H,

% where to store CSD
if ~exist('csdFolder','var') || isempty(csdFolder)
    csdFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\';
end

if ~exist('rebaseline','var')
    rebaseline = 1;
end

% find files in folder
f = what(interpFolder);

%%
for iPP = 1:length(f.mat)
    
    ppID = f.mat{iPP}(1:(regexp(f.mat{iPP}, '_intp')-1));% up until underscore
    
    if ~exist(fullfile(csdFolder, [ppID '_intp_csd.mat']), 'file')
        disp(ppID);
        
        data = load(fullfile(interpFolder, f.mat{iPP}),'erp','RS','ref'); % load 
        
        % un-reference data
        data.erp = data.erp + data.ref;

        tic; % init timer

        data.erp(1:eeg.nChans,:,:) = CSD(data.erp(1:eeg.nChans,:,:),G,H);% do CSD - takes a while
        data.t = toc; % around 15mins for iPP=1
        disp(data.t/60);
        
        % need to re-reference it? Simon's pipeline suggests don't need to
        % so can remove ref then, as it is not CSD'd
        data = rmfield(data, 'ref');


        % re-baseline - sometimes necessary after CSD apparently
        if rebaseline==1
            data.baseline = nanmean(data.erp(:,isBetween(eeg.epochTimes, eeg.blWin),:),2);
            data.erp = data.erp - data.baseline;
        end

        data.isCSD = 1; % set flag
        
        save(fullfile(csdFolder, [ppID '_intp_csd.mat']),'-v7.3',...
            '-struct','data'); % save
    end
end

