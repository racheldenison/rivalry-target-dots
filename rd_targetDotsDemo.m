% rd_targetDotsDemo.m
%
% Demo of the "target dots" method of measuring partial dominance. Rivalry
% stimuli appear for a variable duration. Then the target dots (local
% contrast decrement targets) appear briefly, and the stimulus is removed.
% Participants' task is to report which target dots were visible at the
% time of stimulus offset. Participants may also report their dominant
% percept continuously for the duration of the rivalry presentation.
%
% Modified from rd_testEyeDominance.
%
% Rachel Denison, 28 Feb 2014

global pixelsPerDegree;
global spaceBetweenMultiplier;

%% Set up
% testing location
location = 'laptop';

% i/o
dataDirectoryPath = 'data/';

% params
spatialFreq = 3;
responseDuration = 10;

contrastParam = [0.1];
numContrastVariations = length(contrastParam);

eyeLeftTiltParam = [1 2]; % in which eye is the left-tilted grating?
numEyeLeftTiltVariations = length(eyeLeftTiltParam);

eyeTargetDotsParam = [1 2]; % in which eye are the target dots presented?
numEyeTargetDotsVariations = length(eyeTargetDotsParam);

% number of repetitions of each condition
numReps = 3;

% catch trials
% ** should we have both rivalry catch trials and target dots catch trials
% (aka target dots in both eyes)?
proportionCatchTrials = 0.1;

% alignment between trials
alignmentOption = 'annulusOnly'; % [smiley annulusOnly]
alignmentTargetKeypress = 0; % 1 to wait for a keypress after each alignment target presentation, 0 to go on automatically
alignmentTargetDuration = 1; % time in secons to leave on alignment targets (if not waiting, can set to empty)

% spacing between rivalry images
spaceBetweenMultiplier = 3; 

% sound
soundOn = 1; % 1 for on, 0 for off
v = 1:2000;
soundvector = 0.25 * sin(2*pi*v/30); % a nice beep at 2kH samp freq 

%% Hello to the experimenter
% can write any reminders to the experimenter here
fprintf('Welcome to the Repetition Target Dots study\n\n');

%% Display setup
% PTB-3 correctly installed and functional? Abort otherwise.
AssertOpenGL;

Screen('Preference', 'VisualDebuglevel', 3); % replaces startup screen with black display
screenNumber = max(Screen('Screens'));

% window = Screen('OpenWindow', screenNumber);
window = Screen('OpenWindow', 0, [], [0 0 800 600]);

% Enable alpha-blending, set it to desired blend equation.
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % if using a transparent mask
 
%%%% come back %%%%
% pixelsPerDegree = ang2pix(x);   %number of pixels in one degree of visual angle
% Screen('LoadNormalizedGammaTable', screenNumber, gammaTable);
 
black = BlackIndex(window);  % Retrieves the CLUT color code for black.
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
gray = (black + white) / 2;  % Computes the CLUT color code for gray.  
    
Screen('TextSize', window, 60);

% Set screen to gray
Screen('FillRect',window, gray);
Screen('Flip', window);

%% Keys
% get device numbers
devNums = findKeyboardDevNumsAtLocation(location);
if isempty(devNums.Keypad)
    error('Could not find Keypad! Please check findKeyboardDevNums.')
end

KbName('UnifyKeyNames');

leftKey = '4'; % this is our default key  
rightKey = '5'; % this is our default key  

leftKeyCode (1:256) = 0;
leftKeyCode (KbName (leftKey)) = 1; % creates the keyCode vector

response = 1;
while ~isempty (response)
    fprintf('The key assigned for LEFT is: %s \n', KbName(leftKeyCode));
    response = input ('Hit "enter" to keep this value or a new key to change it.\n','s');
    if ~isempty (response)
        leftKey = response; 
        leftKeyCode (1:256) = 0;
        leftKeyCode (KbName (leftKey)) = 1; % creates the keyCode vector
    end
end

rightKeyCode (1:256) = 0;
rightKeyCode (KbName (rightKey)) = 1; % creates the keyCode vector

response = 1;
while ~isempty (response)
    fprintf('The key assigned for RIGHT is: %s \n', KbName(rightKeyCode));
    response = input ('Hit "enter" to keep this value or a new key to change it.\n','s');
    if ~isempty (response)
        rightKey = response; 
        rightKeyCode (1:256) = 0;
        rightKeyCode (KbName (rightKey)) = 1; % creates the keyCode vector
    end
end

%% Subject ID and filename
subjectID = input ('\nPlease input the alphanumeric code specifier for the subject: ','s');
timeStamp = datestr(now);
fileName = [dataDirectoryPath, subjectID, '_RivalryTargetDots_', datestr(now,'ddmmmyyyy')];

%% Generate trials
% Generate all the possible test configurations and store the values as
% fields in responseArray
trialNum = 1;
for numRep = 1:numReps
    for numContrast = 1:numContrastVariations
        contrast = contrastParam(numContrast);

            for numEyeLeftTilt = 1:numEyeLeftTiltVariations
                eyeLeftTilt = eyeLeftTiltParam(numEyeLeftTilt);

                for numEyeTargetDots = 1:numEyeTargetDotsVariations
                    eyeTargetDots = eyeTargetDotsParam(numEyeTargetDots);

                    responseArray(trialNum).trialType = 0; % regular trial
                    responseArray(trialNum).contrast = contrast;
                    responseArray(trialNum).eyeLeftTilt = eyeLeftTilt;
                    responseArray(trialNum).eyeTargetDots = eyeTargetDots;

                    trialNum = trialNum + 1;

                end % for {numEyeTargetDots}
            end % for {numEyeLeftTilt}
    end % for {numContrast}
end % for numReps
    
% generate catch trials (each run through this loop generates
% two additional trials)
for numCatchTrials = 1:round(trialNum*proportionCatchTrials/2) 
    % create 2 catch trials, grating tilts L & R
    for gratingOrientation = 1:2 % left or right, both eyes' gratings have same orientation
        
        % random contrasts and eye of target dots
        randomNumContrast = round (rand * (numContrastVariations - 1)) + 1;
        contrast = contrastParam(randomNumContrast);
        
        randomNumEyeTargetDots = round (rand * (numEyeTargetDotsVariations - 1)) + 1;
        eyeTargetDots = eyeTargetDotsParam(randomNumEyeTargetDots);
        
        responseArray(trialNum).trialType = 1; % catch trial
        responseArray(trialNum).contrast = contrast;
        responseArray(trialNum).eyeTargetDots = eyeTargetDots;
        responseArray(trialNum).gratingOrientation = gratingOrientation;
        
        trialNum = trialNum +1;
    end 
end

%% Store the experiment info
expt.responseType = responseType;
expt.numReps = numReps;
expt.spatialFreq = spatialFreq;
expt.responseDuration = responseDuration;
expt.alignmentOption = alignmentOption;
expt.soundOn = soundOn;
expt.timeStamp = timeStamp;

%% Ready for the experiment
totalNumTrials = trialNum - 1 % adjust for final increment above

% create random  vector to scramble presentation order
randomOrderVector = randperm (totalNumTrials);

if soundOn
    sound (soundvector, 8000); % make a beep to say we're ready
end

% Opportunity to align stereoscope
presentAlignmentTargetsWaitOption (window, devNums, 1, []) ;  % present alignment targets and wait for a keypress
DrawFormattedText(window, 'Press any key to begin', 0, 'center', [255 255 255]);
Screen('Flip', window);
KbWait(devNums.Keypad); 

%% Present trials
for j = 1:totalNumTrials
    currentTrial = randomOrderVector(j);
    
    contrast = responseArray(currentTrial).contrast;
    eyeTargetDots = responseArray(currentTrial).eyeTargetDots;
    
    switch responseArray(currentTrial).trialType
        case 0 % normal trial
            eyeLeftTilt = responseArray(currentTrial).eyeLeftTilt;

            if eyeLeftTilt == 1
                leftTargetOrient = 1; % left
                rightTargetOrient = 2; % right
            else
                leftTargetOrient = 2; % right
                rightTargetOrient = 1; % left 
            end
        case 1 % catch trial
            gratingOrientation = responseArray(currentTrial).gratingOrientation;
            leftTargetOrient = targetOrientation;
            rightTargetOrient = targetOrientation;
    end
    
    switch alignmentOption
        case 'smiley'
            presentAlignmentTargetsWaitOption (window, devNums, alignmentTargetKeypress, alignmentTargetDuration) ;  % present alignment targets
        case 'annulusOnly'
            if j==1
                presentBlankAlignmentTargetsWaitOption (window, devNums, alignmentTargetKeypress, alignmentTargetDuration) ;  % present alignment targets
            else
                WaitSecs(alignmentTargetDuration);
            end
    end    
    
    if soundOn
        sound (soundvector, 8000); % make a beep
    end
    
    % Here is the trial presentation
    [responseArray(currentTrial).times, responseArray(currentTrial).keyboardEvent]  = ...
        presentRivalryTargetDots (window, spatialFreq, contrast, ...
        eyeTargetDots, leftTargetOrient, rightTargetOrient, ...
        responseDuration, leftKeyCode, rightKeyCode, devNums);

    % Save after each trial
    save(fileName, 'subjectID', 'expt', 'randomOrderVector', 'responseArray');   
end
     
%% Clean up
Screen('CloseAll');  
ShowCursor;
 