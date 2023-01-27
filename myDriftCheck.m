function myDriftCheck(el, fixWinSize, deg2px, waitTime, center, fixCross, drawRed)
% function myDriftCheck(el, fixWinSize, deg2px, waitTime, center, fixCross)
% checks if mean eye-position over 50 samples is within fixWinSize tolerance 
% from center of screen, returns if true at any point during waitTime.
% otherwise, it stops recording, runs a DriftCorrection via Eyelink, 
% applies the correction, restarts recording, and prints the correction.
%
% Inputs:
%   el = structure from ELsetupCalib or EyelinkInitDefaultsSK, containing:
%           el.window: PTB window pointer
%           el.eye_used: 0=left,1=right, -1=null
%           el.yellow: yellow fix cross colour
%           el.foregroundcolor
%           el.backgroundcolour 
%           el.calibrationtargetsize
%           el.calibrationtargetwidth
%   fixWinSize: size of fixation cross in vis deg
%   deg2px: number of pixels in a degree
%   waitTime: number of seconds to wait for mean gaze to be inside fixcross
%   center: [x y] coords of center of screen
%   fixCross: function handle for fixation cross drawing
%   drawRed: 1=draw red square for current eye-position during wait. 0=don't
% 
% 
% John Grogan, 2022. Based on code from Elaine Corbett and Sanjay Manohar.

if ~exist('drawRed', 'var') || isempty(drawRed)
    drawRed = 0;
end

t = GetSecs;
% wait a bit (3s?)
% getting eye position on each sample
x = []; y = [];
while (GetSecs-t) < waitTime
    if (Eyelink('NewFloatSampleAvailable')==1)
        s = Eyelink('NewestFloatSample');

        x = [x; s.gx(el.eye_used+1) - center(1)];
        y = [y; s.gy(el.eye_used+1) - center(2)];
        if size(x,1)>50 % only keep 50 samples
            x = x(end-49:end);
            y = y(end-49:end);
        end

        if numel(x)>1 && drawRed
            Screen('fillrect', el.window, [255 0 0],[x(end) y(end) x(end) y(end)] + [-2 -2 2 2] + [center center]);
%             Screen('fillrect', el.window, [128 128 0], [center center] + [-5 -5 5 5]);
            fixCross(el.window, el.yellow); % fixation cross
            Screen('flip', el.window);
        end
        if mean(sqrt(x.^2 + y.^2)./deg2px) < fixWinSize % if within fix tol
            return;
        end

    end
    WaitSecs(0.001);%
end

% if here, then not fixating

Eyelink('StopRecording'); % pause recording during this

% make target to put around fixation
[width, ~]=Screen('WindowSize', el.window);
sz=round(el.calibrationtargetsize/100*width);
inset=round(el.calibrationtargetwidth/100*width);

rect1=CenterRectOnPoint([0 0 sz sz], center(1), center(2));
rect2=InsetRect(rect1, inset, inset);
                
status = 27;
while status~=0
    Screen( 'FillOval', el.window, el.foregroundcolour,  rect1 ); % big target
    Screen( 'FillOval', el.window, el.backgroundcolour, rect2 ); % small target
    fixCross(el.window, el.yellow); % fixation cross
    Screen( 'Flip',  el.window);

    status = Eyelink('DriftCorrStart', center(1), center(2), 1, 0, 1);  % get the drift correction
end
% final 3 args: 
% dtype=1, do_drift_corr() (0=eyelink_driftcorr_start())
% dodraw=1, eyelink draws a target
% allow_setup=1, can exit to setup in eyelink

% % display these
% s = Eyelink('NewestFloatSample');
% x2 = s.gx(el.eye_used+1) - center(1);
% y2 = s.gy(el.eye_used+1) - center(2);
% disp([mean(x) mean(y); x2 y2]);

Eyelink('ApplyDriftCorr'); % actually apply the correction

% display correction applied
[~,msg] = Eyelink('CalMessage'); % [offset in degrees, x pix, y pix]
disp(msg);


Eyelink('StartRecording'); % restart recording

WaitSecs(0.1); % small pause

