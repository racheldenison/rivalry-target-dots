function selectedTargets = presentTargetDotResponseScreen(window, isi, nTargets, blankTexture, targetRects, white)

% ISI before putting up the response circles
WaitSecs(isi); % Note this is not timed precisely!

% Make blank annulus tex using previously generated image
blanktex = Screen('MakeTexture', window, blankTexture);

% Show response circles
Screen('DrawTexture', window, blanktex);
Screen('FrameOval', window, white, targetRects', 2);
Screen('Flip', window);

% Collect responses
iClick = 1;
keepGoing = 1;
while keepGoing
    % Collect mouse click
    [mouseX, mouseY, buttons] = GetMouse(window);
    if any(buttons([2 3])) % right or middle click
        keepGoing = 0;
        break;
    end
    if buttons(1) % left click
        clickX = mouseX;
        clickY = mouseY;
    else
        clickX = -1; % impossible values for no click
        clickY = -1;
    end
    
    % Find the selected target
    inRectX = clickX > targetRects(:,1) & clickX < targetRects(:,3);
    inRectY = clickY > targetRects(:,2) & clickY < targetRects(:,4);
    selectedCircle = find(inRectX & inRectY);
    
    if ~isempty(selectedCircle) % might be empty if you accidentally click outside the target rects
        % Account for equivalent left and right targets
        selectedTarget = mod(selectedCircle, nTargets); 
        if selectedTarget==0
            selectedTarget = nTargets;
        end
                
        % Update the selected targets list
        selectedTargets(iClick) = selectedTarget;
        iClick = iClick + 1;
        
        % Fill in the selected circle (on both sides)
        selectedRects = targetRects([selectedTargets selectedTargets+nTargets],:);
        Screen('DrawTexture', window, blanktex);
        Screen('FillOval', window', [0 1 0]*white, selectedRects');
        Screen('FrameOval', window, white, targetRects', 2);
        Screen('Flip', window);
    end
end

% Just in case anyone clicks a target twice ...
selectedTargets = unique(selectedTargets);

% Return to blank annulus
Screen('DrawTexture', window, blanktex);
Screen('Flip', window);
