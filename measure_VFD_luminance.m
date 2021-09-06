%% Initialze PTB (PsychToolBox) screen & miscellaneous graphic options
sca;
bg_val = 0; 
Screen('Preference', 'SkipSyncTests', 1); % skip PTB sync-test
Screen('Preference', 'SkipSyncTests', 2);
Screen('Preference', 'VisualDebuglevel', 3); % skip PTB warning
screenNumber = max(Screen('Screens'));  % get the pointer to the screen (select secondary screen in case of multiple displays)
white = WhiteIndex(screenNumber); % define white (usually 255)
grey=GrayIndex(screenNumber); % define gray (usually 127)
[w, windowRect] = Screen('OpenWindow', screenNumber, bg_val); % open the window object and get a pointer to access that window
Screen('Flip', w); % initial flip to clean the screen
ifi = Screen('GetFlipInterval', w); % get interframe interval
hz = 1/ifi; % get refresh rate
MaxPriority(w); % set top priority
[xCenter, yCenter] = RectCenter(windowRect); % get central coordinates (resolution-dependent)
[xScreen, yScreen]= Screen('WindowSize', w); % get screen coordinates (resolution-dependent)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % enable alpha blending to smooth the edge of the stimulus and the VFD conditions
Screen('TextSize', w, 30); % set default text size to 30
AssertOpenGL; % ensure that OpenGL is running properly before running the data acquisition
%%

% create a circle in the center, control its value by clicking
% create a staircase VFD and move it with the mouse to be overlapped with
% the circle in the center

circle_value = 0;
circle_rect = [0 0 100 100];
circle_rect = CenterRectOnPoint(circle_rect, xCenter, yCenter);

mask_value = 0;
maskTest = ones(200,200)*bg_val;
maskTest(:,:,2) = zeros(200,200);
texTest = Screen('MakeTexture', w, maskTest);
dstRectTest=Screen('Rect', texTest);
dstRectTest = CenterRectOnPoint(dstRectTest, xCenter, yCenter);

while 1
    [x,y,buttons] = GetMouse; WaitSecs(0.1);
    Screen('FillRect', w, bg_val, windowRect);
    Screen('FillOval', w, circle_value, circle_rect);
    Screen('DrawText', w, ['circle value: ' num2str(circle_value) ' | mask value: ' num2str(mask_value)], xScreen - 500, yScreen-50, [255 255 255]);
    Screen('DrawTexture', w, texTest,[],dstRectTest);
    if buttons == [1 0 0 0 0]
        if circle_value < 250; circle_value = circle_value+10; end
    elseif buttons == [0 1 0 0 0]
        if circle_value > 0;circle_value = circle_value-10; end
    elseif buttons == [0 0 1 0 0]
        sca; break;
    elseif buttons == [0 0 0 1 0]
        % decrease mask
         if mask_value > 0; mask_value = mask_value - 10; end
        maskTest(:,:,2) = zeros(200,200)+mask_value;
        texTest = Screen('MakeTexture', w, maskTest);
    elseif buttons == [0 0 0 0 1]
        % increase mask
        if mask_value < 250; mask_value = mask_value + 10; end
        maskTest(:,:,2) = zeros(200,200)+mask_value;
        texTest = Screen('MakeTexture', w, maskTest);
    end
   Screen('Flip',w);
end