% Script to convert HFA data into PTB textures to be used as gaze-contingent simulated
% visual field defects.
% It requires: 
% (1) HFA data formatted as a 52x1 array representing the sensitivity values of the visual field
% (2) screen luminance measurement matrix for all combinations of Psychtoolbox (PTB) values vs. mask transparency values (from 0 to 250 in steps of 10)
% (3) physical dimensions of the display screen and viewing distance
%
% ---DEFINITIONS---
% Luminance: the physical luminance measured with a photometer and
% expressed in cd/m2 or asb (this script uses asb, but you can convert
% between the two since 1 cd/m2 = 3.1416 asb)
%
% Value: the scalar value used for PsychtoolBox paramters
%
% Intensity: the perceived luminance expressed in dB. This value depends on
% the general characteristic of the setup on which the measurements are
% made. The stimulus intensity is defined as:
%
% dB = 10 * log(max_lum_asb / (stim_asb - bg_asb))
%
% where 'max_lum_asb' is the maximum luminance value in asb allowed by the
% setup, 'stim_asb' is the luminance in asb of the stimulus and 'bg_asb' is the
% luminance in asb of the background. The same equation can be used to
% calculate subject sensitivity in dB by substituting the stimulus luminance with
% the stimulus at threshold luminance

%% Clean workspace and import required data
clear all;
close all;
clc; disp('Import data...');
arch = readtable('/Volumes/GoogleDrive/Shared drives/R&D/Product Development/scotoma_simulations/archetype_supplementary_materials/supplement/archetypes.xls'); % scotoma archetypes values
HFA = readtable('/Volumes/GoogleDrive/Shared drives/R&D/Product Development/scotoma_simulations/archetype_supplementary_materials/supplement/totaldeviations.xls'); % HFA database values
arch = arch{:,:}; % vectorize the table
HFA = HFA{:,:}; % vectorize the table
screen_lum = readtable('/Volumes/GoogleDrive/Shared drives/R&D/Product Development/scotoma_simulations/screen_luminance_measurement.xlsx','Range','AB4:BA29','ReadVariableNames', false,'Sheet',2); % read specified range of the file containing the luminance matrix in apostillb (asb)
screen_lum = screen_lum{:,:}; % vectorize the table
ptb_val = 0:10:250; % define the ranges over which the measurements were taken
mask_val = 0:10:250;
monit_width = 540; % set the monitor physical width in millimiters
view_dist = 600; % set the viewing distance in millimiters
clc; disp('Import data...DONE');
%% Calculate background and stimulus parameters
% Construct the function F1 to map luminance values onto
% PTB value. We take only the first column of screen
% luminance measurements where mask = 0;
clc; disp('Calculate stimulus and background parameters...');
F1 = csapi(screen_lum(:,1),ptb_val);
bg_asb = 31.5; % set the background to 31.5 asb (standard in many perimeters)
max_lum_asb = 1000; % set the maximum luminance to 1000 asb (some monitors can go higher, but this should be a safer bet to have standardized stimuli)
min_stim_asb = max_lum_asb / exp(32.8/10) + bg_asb; % inverse intensity formula to calculate the minimum stimulus luminance in asb correspondent to foveal sensitivity (approx. 32.8 dB)
max_stim_asb = max_lum_asb; % the maximum stimulus luminance is the same as the maximum allowed by our setup 
bg_val = fnval(F1,bg_asb); % calculate the background PTB value corresponding to desired luminance
max_val = fnval(F1,max_lum_asb); % calculate the maximum available PTB value corresponding to desired luminance (this will be used for the highest intensity stimulus)
stim_val = max_val; % set the stimulus PTB value as the maximum available in our setup
min_stim_val = fnval(F1, min_stim_asb); % calculate the minimum stimulus PTB value for foveal sensitivity
clc; disp('Calculate stimulus and background parameters...DONE');
%% Initialze PTB screen & miscellaneous graphic options. 
clc; disp('Initialize PTB window...');
Screen('Preference', 'SkipSyncTests', 1); % skip PTB sync-test. We do not need ultra-high temporal precision and it might fail on most common screens
Screen('Preference', 'SkipSyncTests', 2);
Screen('Preference', 'VisualDebuglevel', 3); % skip PTB warning
screenNumber = max(Screen('Screens'));  % get the pointer to the screen (select secondary screen in case of multiple displays)
white = WhiteIndex(screenNumber); % define white (usually 255)
grey=GrayIndex(screenNumber); % define gray (usually 127)
[w, windowRect] = Screen('OpenWindow', screenNumber, bg_val); % open the window object and get a pointer to access that window, initialize the window with a background corresponding to 31.5 asb
Screen('Flip', w); % initial flip to clean the screen
ifi = Screen('GetFlipInterval', w); % get interframe interval
hz = 1/ifi; % get refresh rate
MaxPriority(w); % set top priority
[xCenter, yCenter] = RectCenter(windowRect); % get central coordinates (resolution-dependent)
[xScreen, yScreen]= Screen('WindowSize', w); % get screen coordinates (resolution-dependent)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % enable alpha blending to allow PTB textures with transparency masks
Screen('TextSize', w, 30); % set default text size to 30
AssertOpenGL; % ensure that OpenGL is running properly. It is needed for the display of simulated visual field defects
clc; disp('Initialize PTB window...DONE');
%% Spatial conversion from HFA to simulated VFD
clc; disp('Spatial conversion HFA into sVFD...');
n = 456; % index number to read the values of a specific archetype (1 to 17) or of an HFA map from the database (1 to 13231) 
% vf = padding_vf(arch(n,:)); to read archetype instead of HFA, uncomment this line and comment the line below
vf=padding_vf(HFA(n,:)); % make the archetype or HFA map a 8x9 matrix by filling-in with NaNs the locations not measured by the HFA
% substitute the filled-in NaNs with neighbouring information computed by an autoregressive model. 
% In the 'fillgaps' function 'MAXLEN' is set to 3 to avoid taking information from non-neighbouring visual field locations and 'ORDER' is set to 1 to avoid model overfitting. 
% The autoregressive model is applied to the transpose since it prioritizes columns instead of rows (we want the opposite, given the orientation of the original HFA map)
vffilled = fillgaps(vf',3,1); 
vffilled = vffilled'; % re-transpose to put the map back in original format
upvf = imresize(vffilled, [yScreen xScreen], 'bicubic'); % upscale the vf map to match the screen resolution. Now we know that the horizontal resolution corresponds to 9 slots of 5 degrees each (the original size of the HFA map)
% Calculation of scaling factors. They need to be calculated for horizontal and
% vertical axes separatedly since the aspect ratio of the original HFA map
% is different from the aspect ratio of a computer screen
hfa_slot_h = size(upvf,2)/size(vf,2); % ratio between the screen resolution (horizontal) and number of horizontal locations of the HFA. This gives how many pixels correspond to the original 5 degrees of each HFA location
hfa_slot_deg_h = pix2deg(hfa_slot_h,xScreen,monit_width,view_dist); % size in degrees of corresponding to an hfa_slot (without scaling yet), at a viewing distance set by the user 
h_scaling_factor = 5/hfa_slot_deg_h; % % scaling factor to apply to the horizontal screen resolution so that at 600mm of viewing distance an hfa_slot on video corresponds to 5 degrees (the original size in the HFA)
% apply the same calculations to the vertical axis
hfa_slot_v = size(upvf,1)/size(vf,1); 
hfa_slot_deg_v = pix2deg(hfa_slot_v,xScreen,monit_width,view_dist);
v_scaling_factor = 5/hfa_slot_deg_v;
scaled_svfd = imresize(upvf,[size(upvf,1)*v_scaling_factor size(upvf,2)*h_scaling_factor]); % apply the scaling factors to vertical and horizontal dimensions separately. Now the aspect ratio of the sVFD matches that of an HFA, and at a viewing distance set by the user it has the same dimension in degrees
% since the mask will be controlled in a gaze-contingent way, we need to
% enlarge it, without disrupting the scaling done so far, in such a way
% that the subject will be free to look around in every corner of the
% screen. 
h_limit = 2*xScreen; % we set the new dimensions as twice the original mask dimensions
v_limit = 2*yScreen;
% add NaN values on the top and on the sides, in such a way that the new
% mask will be twice as big as the original. 
% NB: this is NOT a rescaling, the aspect ratio of the centermost part stays the same. This is to ensure
% that the spatial property of the simulated VFD are still in line with those of the HFA
top_nanpad = nan(round((v_limit - size(scaled_svfd,1))/2),size(scaled_svfd,2));% TOP padding
side_nanpad = nan(size(scaled_svfd,1)+2*size(top_nanpad,1), round((h_limit - size(scaled_svfd,2))/2));% SIDE padding
svfd = [top_nanpad; scaled_svfd; top_nanpad]; svfd = [side_nanpad svfd side_nanpad]; % apply the padding
svfd=fillgaps(svfd); % analogously as before, we apply another autoregressive model to fill-in the extra space, both horizontally and vertically.  
svfd=fillgaps(svfd');
svfd=svfd'; % re-transpose to put back in original format. Now we have an appropriately scaled mask, covering twice the area of the display monitor, still containing the sensitivity values in dB that need to be converted into PTB values
clc; disp('Spatial conversion HFA into sVFD...DONE');
%% Sensitivity conversion from HFA (dB) to PTB values
clc; disp('Sensitivity conversion HFA into sVFD...');
% calculate the stimulus intensity in dB 
stim_intensity = 10 * log(max_lum_asb / (max_stim_asb - bg_asb));
% plug this value into the equation: L_svfd = L_max / exp((S_hfa + S_stim)/10) + L_bg 
% to convert the dB value of the sVFD map into luminance (asb), taking into
% account the perceived stimulus intensity and background level
lum_svfd_asb = max_lum_asb ./ exp((abs(svfd) + stim_intensity)/10) + bg_asb;
% Construct a 2-D griddedInterpolant F2 from gridded data X,Y,V to map the
% function F1 that links together PTB value (X), mask transparency (Y) and
% luminance (V). This will be used to convert HFA sensitivities (dB) into
% on-screen values for PTB masks
[X,Y] = ndgrid(ptb_val,mask_val); % build a grid from the PTB and transparency mask values used during the screen luminance measurement
V = screen_lum;
F2 = griddedInterpolant(X,Y,V); 
[Xq_gridded,Yq_gridded] = ndgrid(stim_val,mask_val);% Evaluate F2 at gridded query points defined by the stimulus value (Xq) and all the mask values (Yq)
Vq = F2(Xq_gridded,Yq_gridded); % calculate the luminance values corresponding to a stimulus of a given intensity that is overlapped with a PTB mask
F3 = csapi(Vq,mask_val); % interpolate between stimulus luminance and transparency mask values
final_svfd = fnval(F3,lum_svfd_asb); % calculate the final PTB transparency mask values resulting in the desired luminance for a stimulus of a given intensity
clc; disp('Sensitivity conversion HFA into sVFD...');
%% Plot luminance estimation (uncomment section)
% mesh(X,Y,V); hold on; % Plot the gridded data X,Y,V and the interpolated values Vq
% plot3(Xq_gridded(:),Yq_gridded(:),Vq(:),'r.','markersize',20)
% legend('Measured Luminance','Stimulus Luminance','location','best')
% xlabel('Stimulus Value (PTB)'); 
% ylabel('Transparency Mask Value (PTB)');
% zlabel('Combined Luminance (asb)');
%% Build PTB texture
clc; disp('Build PTB texture...');
mask_svfd = ones(size(final_svfd))*bg_val; %bg_val; % set the background layer of the mask equal to the background value
mask_svfd(:,:,2) = final_svfd; % set the transparency layer with the values calculated earlier
texVFD=Screen('MakeTexture', w, mask_svfd); % create a PTB texture
dstRectVFD=Screen('Rect', texVFD); % create a rect associated with the PTB texture
stim_rect = [0 0 50 50]; 
stim_rect = CenterRectOnPoint(stim_rect, xCenter, yCenter);
clc; disp('Build PTB texture...DONE');
%% Display on screen
SetMouse(xCenter,yCenter,w);
while 1
    [x_mouse,y_mouse,buttons] = GetMouse(w); 
    Screen('FillRect', w, bg_val, windowRect);
    Screen('FillOval', w, stim_val, stim_rect);
    dstRectVFD=CenterRectOnPoint(dstRectVFD,x_mouse,y_mouse); % the position of the VFD will be controlled by mouse
    Screen('DrawTexture', w, texVFD,[],dstRectVFD); % draw the VFD on the screen
    if buttons(1) == 1
         sca; break;
    end
   Screen('Flip',w);
   clc; disp('Running.'); 
   clc; disp('Running..'); 
   clc; disp('Running...');
end



























