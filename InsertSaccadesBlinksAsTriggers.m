function [EEG, eyeTriggers] = InsertSaccadesBlinksAsTriggers(EEG, elFileName, dataFolder, minSaccSizePix, startTrigger)
% function [EEG, eyeTriggers] = InsertSaccadesBlinksAsTriggers(EEG, elFileName, dataFolder, minSaccSizePix, startTrigger)
% detect saccades and blinks using snipSaccades.m from matlib
% (github.com/sgmanohar/matlib) and insert triggers into the EEG.event
% structure to show the timings of these detected eye-movements/blinks.
% Inserts trigger codes at start/end of saccades/blinks.
% 
% Requires: matlib files: snipSaccades.m, findregions.m, sq.m, apply.m,
%               readEDFASC.m. 
%           EYE-EEG package: addevents.m
%           cellRegexpi.m (John Grogan)
% 
% Inputs:
%   EEG: eeglab structure containing EEG.srate
%   elFileName: string: eye-link edf file name stub (no extension)
%   dataFolder: path to folder containing elFileName
%   minSaccSizePix: scalar, minimum saccade size in pixels to include.
%       default is zero (i.e. all saccades) but can set to pixelsPerDegree
%       to only include 1 vis deg saccades
%   startTrigger: string: EEG.event.type trigger code for the start-expt
%       trigger that is used to synchronise the EEG and Eye-tracker
%       triggers.
% 
% Outputs:
%   EEG: eeglab structure with start/end-points of saccades/blinks added as
%       events to the EEG.event and EEG.urevent fields
%   eyeTriggers: cell array containing the trigger codes used for these events
% 
% John Grogan, 2022
% 


%% check inputs

if ~exist('minSaccSizePix','var') || isempty(minSaccSizePix)
    minSaccSizePix = 0; % include all
%     minSaccSizePix = 26.4773; % 1 vis deg @ 60cm dist, 1024px width, 40.5cm width
end

if ~exist('startTrigger','var') || isempty(startTrigger)
    startTrigger = '36'; % from PostPulse task
elseif isnumeric(startTrigger)
    startTrigger = num2str(startTrigger);
    warning('converting startTrigger from number to string');
end

%% define codes

eyeTriggers = {'saccadeStart', 54;
               'saccadeEnd',   55;
               'blinkStart',   56;
               'blinkEnd',     57;  }; % will regexp the names, to deal with l/r

%% read asc file into matlab

if ~exist(fullfile(dataFolder, [elFileName '_edfread.mat']), 'file') % if not already exist, make it

%     r=readEDFASC(fullfile(dataFolder, [elFileName '.edf']), 1, 1); % read in asc file with blink/sacc/fix noted
    
    disp('running quicker version of readEDFASC.m, but it can differ slightly due to how Edf2Mat.m and edf2asc.exe process');
    r = QuickReadEDFASC(fullfile(dataFolder, [elFileName '.edf']));
    save(fullfile(dataFolder, [elFileName '_edfread.mat']),'-struct','r') % save it to avoid re-running 

else % otherwise load it
    r = load(fullfile(dataFolder, [elFileName '_edfread.mat']));
end


%% snip saccades to detect all saccades and blinks

% don't interpolate eye-trace, as now NaN show blinks
[raw, info]=  snipSaccades(r, 'TRIGTRSTART_t', 'TRIGTREND_t', 'minsize', minSaccSizePix, 'interpolate',0);

%% get times of saccades in ET timing

saccStartTimes = info.sRT / 1000; % make into ms from 1st trigger
saccDurations = info.sDur/1000;
saccEndTimes = saccStartTimes + saccDurations;

%% find blink timings from trace
% info.sBlink is only blinks within a saccade, so instead find the
% NaN-stretches within the raw eye-trace

% find blinks from non-interpolated traces
blinkInds = apply(2, @(x) {findregions(x)}, isnan(raw));
blinkInds = permute(nancat(blinkInds),[3,2,1]);

% remove trailing nans
blinkInds(repmat(blinkInds(:,2,:) >= size(raw,2),1,2)) = NaN;

blinkStartTimes = sq(blinkInds(:,1,:)); % take start indices (ms from ITI trigger)
blinkEndTimes = sq(blinkInds(:,2,:)); % end indices

% check that there are the same number of each
if any(isnan(blinkStartTimes) & ~isnan(blinkEndTimes), 'all')
    disp('blink starts and ends do not match up')
    keyboard;
end

blinkDurations = blinkEndTimes - blinkStartTimes; % get duration (ms)

% maybe threshold these blink durations?
minBlinkDur = 50; % ms

validBlinks = blinkDurations >= minBlinkDur; % only keep ones that reach this
blinkStartTimes(~validBlinks) = []; % exclude others
blinkEndTimes(~validBlinks) = [];
blinkDurations(~validBlinks) = [];

% convert to seconds, and add to itiTimes (i.e. seconds from start of expt)
blinkStartTimes = blinkStartTimes/1000;
blinkEndTimes = blinkEndTimes/1000;
blinkDurations = blinkDurations/1000;


%% remove any saccades occuring inside a blink

isInsideBlink = any(saccStartTimes >= blinkStartTimes & saccEndTimes<=blinkEndTimes,1); % entire saccade is inside

% % should I also remove those that start or end inside a blink?
% % no, leave for now, as basically acts to extend blink timing for
% flagging
% isInsideBlink = any(saccStartTimes >= blinkStartTimes & saccStartTimes<=blinkEndTimes); % starts inside
% isInsideBlink = any(saccEndTimes >= blinkStartTimes & saccEndTimes<=blinkEndTimes); % ends inside


saccStartTimes(isInsideBlink) = []; % remove
saccEndTimes(isInsideBlink) = [];
saccDurations(isInsideBlink) = [];


%% need to convert all times into samples for EEG.srate

eegOffset = EEG.event(strcmpi({EEG.event.type}, startTrigger)).latency; % # samples to START trigger


saccStartTimes = round(saccStartTimes * EEG.srate + eegOffset); % already in ms. so convert to samples, and shift
saccEndTimes = round(saccEndTimes * EEG.srate + eegOffset);
saccDurations = round(saccDurations * EEG.srate); % no offsets for duration

% blinks - need to add offset still as that's when 'raw' starts
blinkStartTimes = round(blinkStartTimes * EEG.srate + eegOffset); % already in ms. so convert to samples, and shift
blinkEndTimes = round(blinkEndTimes * EEG.srate + eegOffset); 
blinkDurations = round(blinkDurations * EEG.srate); 




%% using eye-eeg function. puts into event &urevent, re-sorts

% saccades - make matrix of metrics to include
saccMat =  [saccStartTimes', saccDurations', saccEndTimes', ...
    repmat(eyeTriggers{1,2}, length(saccStartTimes),1), info.sAmpl(~isInsideBlink)',...
    info.sSpd(~isInsideBlink)', angle(info.sVec(~isInsideBlink))',...
    real(info.sEndpt(~isInsideBlink))', imag(info.sEndpt(~isInsideBlink))'];

EEG = addevents(EEG, saccMat,...
    {'latency','duration','endtime','edftype','sac_amplitude', 'sac_vmax','sac_angle', 'sac_endpos_x', 'sac_endpos_y'},...
    eyeTriggers{1,1});

% also do end-times of saccades
saccMat(:,4) = eyeTriggers{2,2}; % change trigger code to endSacc
% replace startTimes with endTimes in saccMat column order
EEG = addevents(EEG, saccMat(:,[3 2 3 4:end]),...
    {'latency','duration','endtime','edftype','sac_amplitude', 'sac_vmax','sac_angle', 'sac_endpos_x', 'sac_endpos_y'},...
    eyeTriggers{2,1});

% blinks
EEG = addevents(EEG, [blinkStartTimes, blinkDurations, blinkEndTimes repmat(eyeTriggers{3,2}, length(blinkStartTimes),1)],...
    {'latency','duration','endtime','edftype'}, eyeTriggers{3,1});

% also do the end-times of blinks
EEG = addevents(EEG, [blinkEndTimes, blinkDurations, blinkEndTimes repmat(eyeTriggers{4,2}, length(blinkStartTimes),1)],...
    {'latency','duration','endtime','edftype'}, eyeTriggers{4,1});



end
