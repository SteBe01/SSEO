function EarthSphere1(Rt)

% The function gives the Earth's texture loaded in a plot.
% 
% EarthSphere(Rt)
%
% Input arguments:
% ------------------------------------------------------------------------
%   Rt          [1x1]       Earth mean radius       [km]
%
% Output arguments:
% ------------------------------------------------------------------------
%   []          [figure]    Figure open with the Earth picture loaded
%
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.


% Set the default value for the Earth radius in case of no inputs:
if nargin < 1
 Rt = 6378;
end

% Link to the earth image:
Earth_Image = 'Earth.jpg';

% Background colour:
background_plot = 'w';

hold on;
grid on;
axis equal;

% Axes labels
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% Initial view
view(120,30);

%% Create the Earth surface as a wireframe:

% Number of panels necessary to model the sphere:
N_panels = 180; 

% 3D meshgrid of the sphere points through the 'ellipsoid' function:
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, N_panels);

% Create the globe with the 'surf' function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none','handlevisibility','off');


%% Texturemap the globe:

% Load the Earth image for texturemap:
cdata = imread(Earth_Image);


% Transparency of the globe: 
% if alpha = 1: opaque; if alpha = 0: invisible
alpha = 1;

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

hold on
grid on

% I,J,K axis in cartesian state:
% quiver3(0,0,0,15000,0,0,'k','linewidth', 1, 'handlevisibility', 'off')
% quiver3(0,0,0,0,15000,0,'k', 'linewidth', 1, 'handlevisibility', 'off')
% quiver3(0,0,0,0,0,15000,'k', 'linewidth', 1 , 'handlevisibility', 'off')

end