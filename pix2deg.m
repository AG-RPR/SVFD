% Convert pixels in degrees of visual angle
%
% deg = pix2deg(pix, xres, width, view_dist)
%
% pix:          the quantity in pixels
% xres:         horizontal resolution of the screen
% width:        physical width of the screen
% view_dist:    viewing distance

function deg = pix2deg(pix, xres, width, view_dist)
pix(pix<0) = 0;

deg = rad2deg( atan( (pix/2) / (view_dist*xres/width) ) )*2;

end