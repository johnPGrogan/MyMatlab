function s = QuickReadEDFASC(fname)
% function QuickReadEDFASC(fname)
% shortcut version of readEDFASC that will get r.pos, r.saccade, and the
% trigger names, to allow use in snipSaccades.m (as readEDFASC.m takes ages
% to run).
% uses Edf2Mat.m, cellRegexpi.m
% 
% Inputs:
%   fname is edf file name (incl path)
% 
% Outputs:
%   s is a structure containing (coords are in pixels, samples are number of samples from the first)
%       pos, [nSamples x 4] matrix with [sample xCoord yCoord pupilSize]
% 
%       saccade, [nSaccades x 9] matrix with columns [startSample,
%       endSample, duration, startXCoord, startYCoord, endXCoord,
%       endYCoord, amplitude(degrees), peak velocity (deg/sec)]
% 
%       and then 1 field per trigger (with trial/trig numbers removed)
%       containing the sample at which that trigger last occurred in the
%       expt.
% 
% John Grogan, 2022
% 

 
%% read in asc file

m = Edf2Mat(fname); % struct


%% reformat to meet snipSaccades

% snipSaccades needs s(i).pos, s(i).saccade, and s(i).MESSAGES where
% MESSAGES are the trigger names (with numbers removed from them)
% and (i) is trial (which is 1 for this as we treat the entire block as one
% trial)

s = struct();

% s.pos is [time x y pupil], 1 row per sample

s.pos(:,1) = m.Samples.time - m.Samples.time(1); % samples from start of recording
s.pos(:,2) = m.Samples.posX; % x-pos of eyes, in pix
s.pos(:,3) = m.Samples.posY; % y-pos of eyes, in pix
s.pos(:,4) = m.Samples.pupilSize; % pupil size, A.U.


%% saccade structure

% s.saccade has 1 row per detected saccade (detected by eye-tracker online)
% columns are [startT? endT? duration, xStart, yStart, xEnd, yEnd, ampl(deg), peak vel (d/s)]
% mean conversion of the euclidean dist from x,y gives 26.8805 ppd
% (SD=38.02, range = 25.7451-30.1040).
% not sure why it isn't the same ratio for all?


s.saccade(:,1) = m.Events.Ssacc.time' - m.Samples.time(1);
s.saccade(:,2) = m.Events.Esacc.end'  - m.Samples.time(1);
s.saccade(:,3) = m.Events.Esacc.duration' + 1;
s.saccade(:,4) = m.Events.Esacc.posX';
s.saccade(:,5) = m.Events.Esacc.posY';
s.saccade(:,6) = m.Events.Esacc.posXend';
s.saccade(:,7) = m.Events.Esacc.posYend';
s.saccade(:,8) = round(m.Events.Esacc.hypot,2)';
s.saccade(:,9) = round(m.Events.Esacc.pvel)';


%% messages - times of last occurances of the triggers, with numbers removed


sTrigs = m.Events.Messages.info';

% strip the numbers out
sTrigsStrip = regexp(sTrigs, '(\d+)','split');
sTrigsStrip = cellfun(@(x) strcat(x{:}), sTrigsStrip,'UniformOutput',0);

% keep only trigger ones
isTrig = cellRegexpi(sTrigs, 'TRIG')>0; 

sTrigs1 = unique(sTrigsStrip(isTrig)); % unqiues
sTrigs2 = strcat(sTrigs1, '_t'); % append _t


% find times of last matches
for i = 1:length(sTrigs1)
    lastMatch = find(strcmp(sTrigsStrip, sTrigs1{i}),1,'last');
    s.(sTrigs2{i}) = m.Events.Messages.time(lastMatch) - m.Samples.time(1);
end


end
