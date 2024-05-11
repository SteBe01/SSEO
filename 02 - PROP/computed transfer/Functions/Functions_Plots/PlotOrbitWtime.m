function [th0, thf] = PlotOrbitWtime(rr0, vv0, TOF, mu, color, Title, OrbitName, R, J, cD, AM, W)
%PLOTORBITWTIME Summary of this function goes here
%This function takes as inputs the initial conditions (position and
%velocity vectors), time, gravitational constant and (optionally) the data corresponding to the perturbations
%and gives as output the initial and final true
%anomalies, plus a plot of the corresponding branch of orbit.
%
%INPUT: 
% rr0       : starting position vector      [3X1]         [km]
% vv0       : starting velocity vector      [3X1]         [km/s]
% TOF       : time of flight of the orbit   [1X1]         [s]
% mu        : gravitational constant        [1X1]         [km^3/s^2]
% color     : desidered color of the orbit                [string]        
% Title     : desidered title of the plot                 [string]
% OrbitName : desidered name of the orbit                 [string]
% IF PERTURBED:
% R         : radius of the planet (for drag)             [1X1]  [km]
% J         : value of the zonal armonic perturbation     [1X1]  [-] 
% cD        : coefficient of drag                         [1X1]  [-] 
% AM        : A/M, spaceship's rate of surface over mass  [1X1]  [m^2/km] 
%
%OUTPUT:
% th0       : starting true anomaly          [1X1]             [rad]
% thf       : final true anomaly             [1X1]             [rad]
%
%
%Title input must be expressed as: 'Title' (string)
%OrbitName input must be expressed as: 'OrbitName' (string)
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

%Plot of the celestial body being the focus of the orbit:

if astroConstants(13)*0.99<mu && mu<astroConstants(13)*1.01 %astroConstants(13) %mu=mu_earth
    if TOF>0
        tspan = 0 :  100 : TOF ;
    else 
        tspan = 0 : -100 : TOF ;
    end

    R_rif=6378;       %Earth radius [km]
    EarthSphere1(R_rif); %plot of the Earth
    hold on
elseif astroConstants(4)*0.99<mu && mu<astroConstants(4)*1.01 %astroCostants(4) %mu=mu_sun
    if TOF>0
        tspan = 0 :  12*3600 : TOF ;
    else
        tspan = 0 : -12*3600 : TOF ;
    end
    X=0;
    Y=0;
    Z=0;
    % [x, y, z]=sphere;
    % surf(x*6e7, y*6e7, z*6e7, EdgeColor="y", FaceColor='y', HandleVisibility='off');
    plot3(X, Y, Z, '.', MarkerSize=58, HandleVisibility='off', Color="y"); %, DisplayName='Sun'
    hold on

    R_rif=5*1e7;   %per disegnare i versori dei tre assi e il versore del momento della quantità di moto
end

if size(rr0)==[1, 3]  %the function works with column vectors only!
    rr0=rr0';
end

if size(vv0)==[1, 3]  %the function works with column vectors only!
    vv0=vv0';
end

y0 = [rr0,vv0] ; 
c=length(tspan); %number of points for plot
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );   %tolerances are important to prevent gross solutions
if nargin<8
    [~, rrr] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options) ;  %computing the coordinates of the trajectory (ODE solving)
else
    [~, rrr] = ode113(@(t,y) ode_2bp(t,y,mu,R,J,cD,AM,W), tspan, y0, options) ;  %computing the coordinates of the trajectory (ODE solving)
end

rx=rrr(:,1); %x coordinates of the position vector (column)
ry=rrr(:,2); %y coordinates of the position vector (column)
rz=rrr(:,3); %z coordinates of the position vector (column)


%Computing starting and final true anomalies:
[~, ~, ~, ~, ~, th0]=car2kep(rr0, vv0, mu);
rrf=rrr(c, 1:3);                            %final position vector
vvf=rrr(c, 4:6);                            %final velocity vector
[~, ~, ~, ~, ~, thf]=car2kep(rrf, vvf, mu);

%Computing the starting angular momentum:

h=cross(rr0, vv0);  %starting specific angular momentum: once plotted grants to understand the direction of rotation of the orbit

%Computing the variations of a, e:



%Angular momentum for drawing:

h=h/norm(h);  %normalizzo h in modo che la dimensione non influenzi la scala della rappresentazione 
h=3*R_rif*h;     % "allungo" h in modo che non venga coperto dalla superficie terrestre

hx=[0, h(1)];
hy=[0, h(2)];
hz=[0, h(3)];

%plotting: 

% colors:'hot', 'abyss', 'autumn', 'copper', 'cool', 'sky', 'summer',
% 'bone'

if  strlength(color)>1
    rz(c)=NaN;
    C=linspace(0, 100, c);
    p = patch(rx, ry, rz, C, 'FaceColor','none','EdgeColor', 'interp', 'linewidth', 1, DisplayName=OrbitName); %use patch   %'FaceColor','none','EdgeColor','interp'
    CB=colorbar();                                                        %'Direction', 'reverse'
    CB.Label.String = 'ToF% in the orbit';

    oldcmap = colormap(color);                                            %colormap(p, hot); %autumn(100)
    colormap( flipud(oldcmap) );
    hold on
    plot3(hx, hy, hz, 'k.-', HandleVisibility='off'); %Angular momentum DisplayName='Specific Angular Momentum'
else
    width = 1;
    plot3(rx, ry, rz, '-', Color=color, DisplayName=OrbitName, LineWidth = width); %trajectory
    hold on
    %plot3(hx, hy, hz, '.-', Color=color, HandleVisibility='off'); %Angular momentum DisplayName='Specific Angular Momentum'
end

quiver3(0,0,0,2*R_rif,0,0,'k','linewidth', 1, 'handlevisibility', 'off');   %versore I
quiver3(0,0,0,0,2*R_rif,0,'k', 'linewidth', 1, 'handlevisibility', 'off');  %versore J
quiver3(0,0,0,0,0,2*R_rif,'k', 'linewidth', 1 , 'handlevisibility', 'off'); %versore K

%adding Marker for starting & final point:

plot3(rr0(1),rr0(2), rr0(3), '.', 'MarkerSize', 15, MarkerEdgeColor='g', HandleVisibility='off'); %, DisplayName='Starting Point'
plot3(rrf(1),rrf(2), rrf(3), '.', 'MarkerSize', 15, MarkerEdgeColor='r', HandleVisibility='off'); %, DisplayName='Final Point'
legend();

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title(Title);
view(60, 60);
axis equal
axis vis3d
grid on

    

