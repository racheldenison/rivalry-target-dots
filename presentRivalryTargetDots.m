function [responseTimes, responseKeyboardEvent] = ...
    presentRivalryTargetDots(window, spatialFrequency, contrast, ...
        eyeTargetDots, leftTargetOrient, rightTargetOrient, ...
        responseDuration, leftKeyCode, rightKeyCode, devNums)

% presentRivalryTargetDpts
%
% The default targets are 1.8 degrees in diameter, with vergence annulus of
% 2.6 degrees in diameter. This can be easily modified.
%
% output parameters
%   responseArray is an array with a record of times and buttons pushed
%
% Modified from presentRivalryTargetsRepetition_ContinuousResponse.m
%
% Rachel Denison - 28 February 2014


global pixelsPerDegree;

if isempty(pixelsPerDegree)
    pixelsPerDegree = 99; % this is with NEC Monitor with subject sitting 5 feet from the screen. 1280 x 1024 pixel fullscreen.
    display ('in present rivalry targets - passing global variable ppd did not work');
end;

global spaceBetweenMultiplier;
if isempty(spaceBetweenMultiplier)
    spaceBetweenMultiplier = 2;
    display ('error is passing global sbm in present rivalry targs...');
end;

echo off  % Prevents MATLAB from reprinting the source code when the program runs.

pixelsPerCycle = pixelsPerDegree / spatialFrequency; % How many pixels will each period/cycle occupy?
cyclesPerPixel = 1/pixelsPerCycle; % How many periods/cycles are there in a pixel?
radiansPerPixel = cyclesPerPixel * (2 * pi); % = (periods per pixel) * (2 pi radians per period)


% Set diameter of rivalry gratings as well as diameter and thickness of Convergence
% annulus
circleDiameterDegrees = 1.8; % Diameter of presented circle in degrees of visual field.
convergenceAnnulusDiameterDegrees = 2.6; % Diameter of black annulus which individually surrounds both the target gratings
convergenceAnnulusThicknessDegrees = .2;

circleDiameterPixels = circleDiameterDegrees * pixelsPerDegree;
convergenceAnnulusDiameterPixels = convergenceAnnulusDiameterDegrees * pixelsPerDegree;
convergenceAnnulusThicknessPixels = convergenceAnnulusThicknessDegrees * pixelsPerDegree;

widthOfGrid = convergenceAnnulusDiameterPixels ;   %the next lines make sure that it is a whole, even number so matrix indices don't choke.
widthOfGrid = round (widthOfGrid);
if mod (widthOfGrid, 2) ~= 0
    widthOfGrid = widthOfGrid + 1 ;
end;

halfWidthOfGrid =  (widthOfGrid / 2);
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  % widthArray is used in creating the meshgrid.

% Set response target dimensions
responseCircleDiameterDegrees = 0.5;
responseCircleDiameterPixels = responseCircleDiameterDegrees * pixelsPerDegree;

% Set tilt (angle of grating) for Target1 and Target2
tiltInDegreesOrient1Target = 135; % The tilt of the grating in degrees. (135 deg = left)
tiltInDegreesOrient2Target = 45; % The tilt of the grating in degrees. (45 deg = right)
tiltInRadiansOrient1Target = tiltInDegreesOrient1Target * pi / 180; % The tilt of the grating in radians.
tiltInRadiansOrient2Target = tiltInDegreesOrient2Target * pi / 180; % The tilt of the grating in radians.

[width, height]=Screen('WindowSize', window); % in case we need this for a future version...
cx = round(width/2);
cy = round(height/2);

% ---------- Color Setup ----------
% Gets color values.

% Retrieves color codes for black and white and gray.
black = BlackIndex(window);  % Retrieves the CLUT color code for black.
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
gray = (black + white) / 2;  % Computes the CLUT color code for gray.

% Taking the absolute value of the difference between white and gray will
% help keep the grating consistent regardless of whether the CLUT color
% code for white is less or greater than the CLUT color code for black.
absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

% ---------- Image Setup ----------
% Stores the image in a two dimensional matrix.

% Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
% the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
% and y = y(x0, y0) corresponds to the y-coordinate of element "i"
[x y] = meshgrid(widthArray, widthArray);

% Now we create a circle mask
circularMaskMatrix = (x.^2 + y.^2) < (circleDiameterPixels/2)^2;

%Now we create an annulus that surrounds our circular grating
convergenceAnnulus = ...
    ((x.^2 + y.^2) >= (convergenceAnnulusDiameterPixels/2 - convergenceAnnulusThicknessPixels)^2 )  & ...
    ((x.^2 + y.^2) <  (convergenceAnnulusDiameterPixels/2)^2 );

% if we want dot in middle: ...|((x.^2 + y.^2) <  (pixelsPerDegree *.05)^2 ) ;

convergenceAnnulus = ~convergenceAnnulus; %when we multiply this, it will create an annulus of zero/black

responseCircle = (x.^2 + y.^2) < (responseCircleDiameterPixels/2)^2;
responseCircle = ~responseCircle;

%----------------------------------------------------------------------
michelsonContrast = contrast;

a1 = cos(tiltInRadiansOrient1Target) * radiansPerPixel;
b1 = sin(tiltInRadiansOrient1Target) * radiansPerPixel;
phase = 0;

imageMatrixOrient1Target =...
    (gray + absoluteDifferenceBetweenWhiteAndGray * michelsonContrast * ...
    sin(a1*x+b1*y+phase).* circularMaskMatrix)  .*  convergenceAnnulus;


a2=cos(tiltInRadiansOrient2Target)*radiansPerPixel;
b2=sin(tiltInRadiansOrient2Target)*radiansPerPixel;

imageMatrixOrient2Target= ...
    (gray + absoluteDifferenceBetweenWhiteAndGray* michelsonContrast * ...
    sin(a2*x+b2*y-phase).* circularMaskMatrix)  .* convergenceAnnulus;

%----------------------------------------------------------------------
% Note that each entry of matrices varies between minus one and one;
% multiplying these matrices by absoluteDifferenceBetweenWhiteAndGray
% and adding the gray CLUT color code baseline
% converts each entry of imageMatrix into a shade of gray:
% minus one is black; zero, is gray; one is white.
graySpacerMatrix =  ones(widthOfGrid+1,(widthOfGrid)*spaceBetweenMultiplier ) * gray;

blankTargetMatrix = ones(widthOfGrid+1, widthOfGrid+1) * gray;
blankedGratingMatrix = blankTargetMatrix .*convergenceAnnulus;

blankTexture = [blankedGratingMatrix, graySpacerMatrix, blankedGratingMatrix];
blanktex = Screen('MakeTexture', window, blankTexture);

% response target 
responseTargetMatrix = blankTargetMatrix .*convergenceAnnulus .*responseCircle;
responseTexture = [responseTargetMatrix, graySpacerMatrix, responseTargetMatrix];
responsetex = Screen('MakeTexture', window, responseTexture);

% Build rivalry grating textures:
if leftTargetOrient == rightTargetOrient % catch trial, identical target gratings
    if leftTargetOrient == 1
        gratingTexture = [imageMatrixOrient1Target, graySpacerMatrix, imageMatrixOrient1Target];
    elseif leftTargetOrient == 2
        gratingTexture = [imageMatrixOrient2Target,graySpacerMatrix, imageMatrixOrient2Target];
    end
else % normal trial, orthogonal target gratings
    if leftTargetOrient == 1
        gratingTexture = [imageMatrixOrient1Target, graySpacerMatrix, imageMatrixOrient2Target];
    elseif leftTargetOrient == 2
        gratingTexture = [imageMatrixOrient2Target, graySpacerMatrix, imageMatrixOrient1Target];
    end
end

tex = Screen('MakeTexture', window, gratingTexture);

% Target dots
% ** should move these out of this function
nTargetDots = 6;
targetColor = 0.5;
targetSz = 50;
targetSigma = 8;
targetAmp = 0.5;

% Find the dot positions
% First find the distance between the center of the screen and the center of the gratings
gratingCenterDistance = (size(imageMatrixOrient1Target,2) + ...
    size(graySpacerMatrix,2))/2;
% Next find the distance between the 
centerTargetDistance = size(imageMatrixOrient1Target,1)*0.5; % put the targets half way between the center of the grating and the edge of the annulus

targetAngle = 2*pi/(nTargetDots-1); % subtract 1, since we will add one in the center
targetAngles = [0:targetAngle:2*pi-targetAngle];
targetXs = centerTargetDistance*sin(targetAngles);
targetYs = centerTargetDistance*cos(targetAngles);
% add a target in the center of the grating
targetXs = [targetXs 0];
targetYs = [targetYs 0];

if eyeTargetDots==1 % target dots in left eye
    targetPositions = [cx - gratingCenterDistance + targetXs; cy + targetYs]';
else    
    targetPositions = [cx + gratingCenterDistance + targetXs; cy + targetYs]';
end

% Make a target dot
% 2-layer masking
targetDot(:,:,1) = ones(targetSz)*targetColor;
targetDot(:,:,2) = make2DGaussianCentered(targetSz, targetSz, 0, 0, targetSigma, targetAmp); % transparent layer, 0 is completely transparent

targettex = Screen('MakeTexture', window, targetDot*white);

% Switch to realtime:
priorityLevel=MaxPriority(window);

% kluge to deal with random intermittent crashes until MacOS is updated
successfullySetPriority = 0;
while ~ successfullySetPriority
    try
        Priority(priorityLevel);
        successfullySetPriority = 1;

    catch
        successfullySetPriority = 0;
    end
end

% Present rivalry sequence:
% Initialize response measures
dataIndex = 1;
responseTimes(1) = 0;
responseKeyboardEvent(1) = 0;
leftKeyWasDown = 0;
rightKeyWasDown = 0;

timeFlipBlank = GetSecs;

% Display the rivalry gratings
% Draw images
Screen('DrawTexture', window, tex);
timeFlipImage = Screen('Flip', window);

% Check for continuous response to final rivalry display
tic % start timer
while GetSecs < timeFlipImage + responseDuration
    WaitSecs(.01);

    %  Check the Subject Responses
    [ keyIsDown, seconds, keyCode ] = KbCheck(devNums.Keypad); % Check the state of the keyboard.

    leftKeyIsNowDown =  isequal(leftKeyCode, (leftKeyCode & keyCode));
    rightKeyIsNowDown =  isequal(rightKeyCode, (rightKeyCode & keyCode));

    if (~leftKeyWasDown) && leftKeyIsNowDown
%        display ('You just PRESSED the LEFT key');
        leftKeyWasDown = 1;
        responseTimes(dataIndex) = toc;
        responseKeyboardEvent(dataIndex) = 1;
        dataIndex = dataIndex +1;
    elseif leftKeyWasDown && ~leftKeyIsNowDown
%        display ('You just RELEASED the LEFT key');
        leftKeyWasDown = 0;
        responseTimes(dataIndex) = toc;
        responseKeyboardEvent(dataIndex) = 2;
        dataIndex = dataIndex +1;
    elseif (~rightKeyWasDown) && rightKeyIsNowDown
%        display ('You just PRESSED the RIGHT key');
        rightKeyWasDown = 1;
        responseTimes(dataIndex) = toc;
        responseKeyboardEvent(dataIndex) = 3;
        dataIndex = dataIndex +1;
    elseif rightKeyWasDown && ~rightKeyIsNowDown
%        display ('You just RELEASED the RIGHT key');
        rightKeyWasDown = 0;
        responseTimes(dataIndex) = toc;
        responseKeyboardEvent(dataIndex) = 4;
        dataIndex = dataIndex +1;
    end %if
end

% mark the end of the response duration in case key is still down when
% the response duration ends
responseTimes(dataIndex) = toc;
responseKeyboardEvent(dataIndex) = 99;
dataIndex = dataIndex + 1;

% make sure we don't leave a lingering image
Screen('DrawTexture', window, blanktex);
Screen('Flip', window);

% Shut down realtime-mode:
% kluge to deal with random intermittent crashes until MacOS is updated
successfullySetPriority = 0;
while ~ successfullySetPriority
    try
        Priority(0);
        successfullySetPriority = 1;

    catch
        successfullySetPriority = 0;
    end
end

% release all textures and offscreen windows
Screen('Close');

end
