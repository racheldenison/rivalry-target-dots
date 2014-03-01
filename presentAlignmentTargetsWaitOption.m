function  presentAlignmentTargetsWaitOption (window, devNums, waitForKeypress, waitTimeInSecs)

% puts identical left and right target images on screen so subjects can 
% adjust the stereoscope to achieve binocular fusion 
% window is pointer to window to present targets in (assumes gray bkgd)
%
% waits for a key press and then returns
%
% by David Egert - 19 Oct 2007

global pixelsPerDegree;

if isempty(pixelsPerDegree)
    pixelsPerDegree = 99; 
    display ('global variable ppd passing did not work in present alignment targs') % this is with NEC Monitor with subject sitting 5 feet from the screen. 1280 x 1024 pixel fullscreen.
end;

global spaceBetweenMultiplier;

if isempty(spaceBetweenMultiplier)
    spaceBetweenMultiplier = 2;  
    display ('error in passing global sbm in present alignment targs...');
end;

% Initializes the program's parameters.

echo off  % Prevents MATLAB from reprinting the source code when the program runs.

% devNums = findKeyboardDevNums; % devNums.Keypad will be what we want...

% Set diameter of Targets as well as diameter and thickness of Convergence
% annulus

circleDiameterDegrees = 1.8; % Diameter of presented circle in degrees of visual field.
convergenceAnnulusDiameterDegrees = 2.6; % Diameter of black annulus which individually surrounds both the target gratings
convergenceAnnulusThicknessDegrees = .2;

circleDiameterPixels = circleDiameterDegrees * pixelsPerDegree;
convergenceAnnulusDiameterPixels = convergenceAnnulusDiameterDegrees * pixelsPerDegree;
convergenceAnnulusThicknessPixels = convergenceAnnulusThicknessDegrees* pixelsPerDegree;

widthOfGrid =  convergenceAnnulusDiameterPixels;
widthOfGrid = round (widthOfGrid);
if mod (widthOfGrid, 2) ~= 0
    widthOfGrid = widthOfGrid + 1;
end;

halfWidthOfGrid =  (widthOfGrid / 2);
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  % widthArray is used in creating the meshgrid.

   
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
        ((x.^2 + y.^2) <  (convergenceAnnulusDiameterPixels/2)^2 )  ;
    convergenceAnnulus = ~convergenceAnnulus; %when we multiply this, it will create an annulus of zero/black
    
    
	%----------------------------------------------------------------------
   
       
    alignmentImageMatrix  = ...   
        ( x == 0 ) | (y == 0) | ... % the cross
        (  (((x - halfWidthOfGrid/3).^2 + (y + halfWidthOfGrid/3  ).^2)  < (halfWidthOfGrid/8  ).^2) ) | ... % one eye
        (  (((x +  halfWidthOfGrid/3).^2 + (y + halfWidthOfGrid/3  ).^2)  < (halfWidthOfGrid/8  ).^2) )  | ...% another eye 
        (( ((x.^2 + y.^2) >= (halfWidthOfGrid/2 - 10 ).^2 )  & ... % this line and the next are the smile...
        ((x.^2 + y.^2) < (halfWidthOfGrid/2 ).^2 )) & (y > 30)); 
                         
    alignmentImageMatrix =...
      (gray + 0.75 * absoluteDifferenceBetweenWhiteAndGray *  alignmentImageMatrix).*  convergenceAnnulus;
      
    graySpacerMatrix =  ones(widthOfGrid+1,(widthOfGrid)*spaceBetweenMultiplier)* gray;
     
    texture =Screen('MakeTexture', window, [alignmentImageMatrix,graySpacerMatrix, alignmentImageMatrix]);     
   
    % Draw gray full-screen rectangle to clear to a defined background color:
    Screen('FillRect',window, gray);
    Screen('Flip', window);% Show the gray background, return timestamp of flip in 'vbl'

    Screen('DrawTexture', window, texture );
    Screen('Flip', window);    % Show alignmentImage on screen 
         
    if waitForKeypress == 1
        KbWait(devNums.Keypad); % wait for a keypress and then return...
    else
        WaitSecs(waitTimeInSecs);
    end
    
 
end
