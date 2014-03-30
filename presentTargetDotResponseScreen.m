function presentTargetDotResponseScreen(window, isi, blankTexture, targetRects, devNums, white)

% ISI before putting up the response circles
WaitSecs(isi); % Note this is not timed precisely!

% Make blank annulus tex using previously generated image
blanktex = Screen('MakeTexture', window, blankTexture);

% Show response circles
Screen('DrawTexture', window, blanktex);
Screen('FrameOval', window, white, targetRects', 2);
Screen('Flip', window);

% Collect responses
KbWait(devNums.Keypad);

% Return to blank annulus
Screen('DrawTexture', window, blanktex);
Screen('Flip', window);
