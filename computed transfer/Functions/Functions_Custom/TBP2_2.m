function [rr_t, vv_t, rr_bi, vvf] = TBP2_2(M, m, rr0, vv0, rr1, vv1, t, dt)
%TBP Summary of this function goes here
% M : mass of the Sun [kg]
% m : mass of the Earth [kg]
% rr0 : initial position of satellite wrt Earth [km]   IN EARTH RTH
% vv0 : initial velocity of satellite wrt Earth [km/s] IN EARTH RTH
% rr1 : initial position of Earth wrt Sun [km]
% vv1 : initial velocity of Earth wrt Sun [km/s]
% t : duration of the mission [days]
% dt : discretization of the mission [s]
% Rl : radius of the Earth

%output:
%
%
% rrf: final position of third body wrt to Earth in Earth rth frame
% vvf: final velocity of third body wrt to Earth in Earth rth frame

G=6.6743e-20; %[km^3/(kg*s)]

t=t*3600*24;
tspan=0:dt:t;
c=length(tspan); %number of points for plot

MU=G*M;
mu=G*m;

%GET rr0, vv0 IN SUN CENTRED INERTIAL FRAME:
r1n=rr1/norm(rr1);
th1=atan2(r1n(2), r1n(1));
ang1=pi-th1;
R1= [cos(ang1), sin(ang1), 0;
    -sin(ang1), cos(ang1), 0;
     0,         0,        1];
rr0=R1*rr0';
vv0=R1*vv0';

%SOLVE THE ODE
y0 = [rr1,vv1] ;         %state of the Earth in sun centred inertial frame
y1 = [rr1+rr0,vv1+vv0] ; %state of the body in sun centred inertial frame

x0 = [y0, y1] ; 

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );   %tolerances are important to prevent gross solutions
[~, xx] = ode113(@(t,x) ode_3bp(t,x,MU,mu), tspan, x0, options) ;  %computing the coordinates of the trajectory (ODE solving)

%PLOT EARTH CENTRED RTH FRAME:
figure(1)
rr_b=xx( :, 7:9)-xx( :, 1:3);   %position of the body wrt to Earth (inertial)
vv_b=xx( :, 10:12)-xx( :, 4:6); %velocity of the body wrt to Earth (inertial)

rE = xx(:, 1:3); %position of Earth wrt Sun
size(rE)
for i=1:c
    rE_n(i, 1:3) = rE(i, 1:3) / norm(rE(i, 1:3)); %versor of the position of Earth wrt Sun

    th_E(i)=atan2(rE_n(i, 2), rE_n(i, 1));  %angle of the position of Earth wrt Sun (th_E=0 for rE_n=[1, 0, 0])

    ang(i)=-pi+th_E(i);

    Ri= [cos(ang(i)), sin(ang(i)), 0;
        -sin(ang(i)), cos(ang(i)), 0;
          0,          0,           1];

    rr_bi(i, :) = Ri* rr_b(i, :)';  %position of body wrt Earth in rth earth centred frame

end

rxi=rr_bi(:, 1);
ryi=rr_bi(:, 2);
rzi=rr_bi(:, 3);

rrf=rr_bi(end, :);
vvf=(Ri* vv_b(end, :)')';

rzi(c)=NaN;
C=linspace(0, 100, c);
plot3(rxi, ryi, rzi, 'color', "#EDB120", DisplayName='SOHO trajectory');
% p = patch(rxi, ryi, rzi, C, 'FaceColor','none','EdgeColor', 'interp', 'linewidth', 1, HandleVisibility='off'); %use patch   %'FaceColor','none','EdgeColor','interp'
% CB=colorbar();                                                        %'Direction', 'reverse'
% CB.Label.String = 'ToF% in the trajectory';
oldcmap = colormap('cool');                                            %colormap(p, hot); %autumn(100)
colormap( flipud(oldcmap) );
hold on
% plot3(rxi(1), ryi(1), rzi(1), '.', 'MarkerSize', 15, MarkerEdgeColor='g', HandleVisibility='off'); %, DisplayName='Starting Point'
hold on
R_rif=6378;       %Earth radius [km]
% EarthSphere1(R_rif); %plot of the Earth
hold on
% Position of L1 at starting:
f = @(b) 1+(mu/MU)/((1-b)*b^2)-1/(1-b)^3 ;
b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)
xLag1=norm(rE(1, 1:3))*(b);%[km] L1 position wrt to Earth
% plot3(xLag1, 0, 0, '.', 'MarkerSize', 15, MarkerEdgeColor='b', DisplayName='L1 approx. point'); % HandleVisibility='off',
hold on

R_rif=1e5;
% quiver3(0,0,0,2*R_rif,0,0,'k','linewidth', 1, 'handlevisibility', 'off');   %versore I
% quiver3(0,0,0,0,2*R_rif,0,'k', 'linewidth', 1, 'handlevisibility', 'off');  %versore J
% quiver3(0,0,0,0,0,2*R_rif,'k', 'linewidth', 1 , 'handlevisibility', 'off'); %versore K
legend();
xlabel('\textbf{Earth-Sun dir. (x) [km]}', Interpreter='latex', FontSize=10)
ylabel('\textbf{y [km]}', Interpreter='latex', FontSize=10)
zlabel('\textbf{h vect. (z) [km]}', Interpreter='latex', FontSize=10)
% title('Earth_r_t_h');
view(60, 60);
axis equal
axis vis3d
grid on

% % plot of the Earth
% X=0;
% Y=0;
% Z=0;
% plot3(X, Y, Z, '.', MarkerSize=58, HandleVisibility='off', Color="g"); %, DisplayName='Sun'
% hold on


%PLOT OF SUN CENTRED INERTIAL FRAME:
figure(2)
rr_t=xx( :, 7:9);      %position of the body wrt to Sun               
vv_t=xx( :, 10:12);    %velocity of the body wrt to Sun            

%

rx=rr_t(:,1); %x coordinates of the position vector (column)
ry=rr_t(:,2); %y coordinates of the position vector (column)
rz=rr_t(:,3); %z coordinates of the position vector (column)

rEx=rE(:,1); %x coordinates of the position vector (column)
rEy=rE(:,2); %y coordinates of the position vector (column)
rEz=rE(:,3); %z coordinates of the position vector (column)

%plot of body wrt Sun
rz(c)=NaN;
C=linspace(0, 100, c);
% p = patch(rx, ry, rz, C, 'FaceColor','none','EdgeColor', 'interp', 'linewidth', 1, DisplayName='SOHO'); %use patch   %'FaceColor','none','EdgeColor','interp'
% CB=colorbar();                                                        %'Direction', 'reverse'
% CB.Label.String = 'ToF% in the orbit';
oldcmap = colormap('hot');                                            %colormap(p, hot); %autumn(100)
colormap( flipud(oldcmap) );
hold on
% plot3(rx(1),ry(1), rz(1), '.', 'MarkerSize', 15, MarkerEdgeColor='g', HandleVisibility='off'); %, DisplayName='Starting Point'
hold on

%plot of Earth wrt Sun
width = 1;
plot3(rEx, rEy, rEz, '-', Color='g', DisplayName='Earth', LineWidth = width); %trajectory
hold on

R_rif=1e8;
% quiver3(0,0,0,2*R_rif,0,0,'k','linewidth', 1, 'handlevisibility', 'off');   %versore I
% quiver3(0,0,0,0,2*R_rif,0,'k', 'linewidth', 1, 'handlevisibility', 'off');  %versore J
% quiver3(0,0,0,0,0,2*R_rif,'k', 'linewidth', 1 , 'handlevisibility', 'off'); %versore K
legend();
xlabel('\textbf{Earth-Sun dir. (x) [km]}', Interpreter='latex', FontSize=10)
ylabel('\textbf{y [km]}', Interpreter='latex', FontSize=10)
zlabel('\textbf{h vect. (z) [km]}', Interpreter='latex', FontSize=10)
title('solar system');
view(60, 60);
axis equal
axis vis3d
grid on

% plot of the Sun
X=0;
Y=0;
Z=0;
    % [x, y, z]=sphere;
    % surf(x*6e7, y*6e7, z*6e7, EdgeColor="y", FaceColor='y', HandleVisibility='off');
plot3(X, Y, Z, '.', MarkerSize=58, HandleVisibility='off', Color="y"); %, DisplayName='Sun'
hold on


