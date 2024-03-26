function plotOrbit(a, e, i, Om, om, th0, thf, mu, p, optn, color, Title, OrbitName, R, J, cD, AM, W)  %metti Tf come input
%PLOTORBIT plotta graficamente l'orbita descritta dai parametri:
% a   : semimajor axis                                            [1x1] [km]
% e   : eccentricity                                              [1x1] [1]
% i   : inclination                                               [1x1] [rad]
% Om  : right ascension of ascending node (RAAN)                  [1x1] [rad]
% om  : argument of pericenter                                    [1x1] [rad]
% th0 : initial true anomaly                                      [1x1] [rad]
% thf : final true anomaly                                        [1x1] [rad]
% dth : anglular discretization for the orbit                     [1x1] [rad]
% mu  : standard gravitational parameter                          [1x1] [km^3/s^2]
% p   : semilatum rectum (for parabola)                           [1x1] [km]
% optn      : if equal to 'n', shows the pericentre and apocentre [string]
% color     : desidered color of the orbit                        [string]        
% Title     : desidered title of the plot                         [string]
% OrbitName : desidered name of the orbit                         [string]
%
%OUTPUT:
% T1      :  time of flight from th0 to thf                       [1x1] [s]
%IF ELLIPTIC ORBIT:
% T2      :  orbital period (from th0=0 to thf=2*pi)              [1x1] [s]

% T=[T1] if e>=1, T=[T1, T2] if e<1

%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.


%the function handles values expressed in the IS  (km, s, etc)

if e==1 && a==Inf

    if nargin<10

        error("Inserire il valore di p = 2*rP negli Input!");

    end
   
    alpha=pi;   %anomalia limite per orbite paraboliche

    alphat=pi/1.05; %anomalia limite "di sicurezza" per plotting
   
    if th0<=-alpha

       th0=-alphat;       % risetto th0 in modo che appartenga all'orbita
       warning("l'anomalia iniziale è inferiore all'angolo limite!");
       angolo_limite=-alpha;
       disp(angolo_limite);
   
    end

    if th0>=alpha

       th0=-alphat;       % risetto th0 in modo che appartenga all'orbita
       warning("l'anomalia iniziale è superiore all'angolo limite!");
       angolo_limite=alpha;
       disp(angolo_limite);
   
    end

    if thf>=alpha

       thf=alphat;       % risetto thf in modo che appartenga all'orbita
       warning("l'anomalia finale è superiore all'angolo limite!");
       angolo_limite=alpha;
       disp(angolo_limite);

   
    end

    if thf<=-alpha

       thf=alphat;       % risetto thf in modo che appartenga all'orbita
       warning("l'anomalia finale è inferiore all'angolo limite!");
       angolo_limite=-alpha;
       disp(angolo_limite);
   
    end

end



if e>1                      %controllo eccentricità

    if a>=0
        error("Il semiasse maggiore di un'orbita iperbolica deve essere negativo in segno!");
    end

    alpha=pi-acos(1/e);   %anomalia limite per orbite iperboliche
    alphat=alpha-0.1;    %anomalia limite "di sicurezza" per plotting 

    if th0<=-alpha
       th0=-alphat;       % risetto th0 in modo che appartenga all'orbita
       warning("l'anomalia iniziale è inferiore all'angolo limite!");
       angolo_limite=-alpha;
       disp(angolo_limite);

    end

    if th0>=alpha
       th0=-alphat;       % risetto th0 in modo che appartenga all'orbita
       warning("l'anomalia iniziale è superiore all'angolo limite!");
       angolo_limite=alpha;
       disp(angolo_limite);

    end

    if thf>=alpha
       thf=alphat;        % risetto thf in modo che appartenga all'orbita
       warning("l'anomalia finale è superiore all'angolo limite!");
       angolo_limite=alpha;
       disp(angolo_limite);

    end

    if thf<=-alpha
       thf=alphat;        % risetto thf in modo che appartenga all'orbita
       warning("l'anomalia finale è inferiore all'angolo limite!");
       angolo_limite=-alpha;
       disp(angolo_limite);

    end
end

if e<1
    rr_max=kep2car(a, e, i, Om, om, pi, mu);   %distanza all'apocentro (se orbita ellittica)
else
    rrmax=max(norm(kep2car(a, e, i, Om, om, th0, mu)), norm(kep2car(a, e, i, Om, om, thf, mu))); %maximum distance in current plot. useful for plotting the asymptots
end

if e~=1  %ellisse o iperbole
    rr_min=kep2car(a, e, i, Om, om,  0, mu);   %distanza al pericentro
    [r,v] = kep2car(a, e, i, Om, om, th0, mu) ; 
    y0 = [r,v] ; 
elseif e<0
    warning("Attenzione, il modulo dell'eccentricità non è mai minore di zero!!" );
elseif e==1 && a==Inf    %parabola
    rr_min=kep2car(a, e, i, Om, om,  0, mu, p);   %distanza al pericentro
    [r,v] = kep2car(a, e, i, Om, om, th0, mu, p) ; 
    y0 = [r,v] ;                
end

%Plot della sfera del corpo celeste attorno cui avviene l'orbita per disegno:

TOF = ToF(a, e, th0, thf, mu);

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

c=length(tspan);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );   % le tolleranze sono importanti se no l'orbita ti viene troppo approssimata
if nargin<14
    [~, rrr] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options) ;  %calcolo delle coordinate sulla traiettoria (risoluzione ODE)
else
    [~, rrr] = ode113(@(t,y) ode_2bp(t,y,mu,R,J,cD,AM,W), tspan, y0, options) ;  %calcolo delle coordinate sulla traiettoria (risoluzione ODE)
end

rx=rrr(:,1); %componenti x del raggio vettore (in colonna)
ry=rrr(:,2); %componenti y del raggio vettore (in colonna)
rz=rrr(:,3); %componenti z del raggio vettore (in colonna)

rr0=rrr(1, :);  %vettore posizione iniziale

rrf=rrr(c, :);  %vettore posizione finale

h=cross(r(:, 1), v(:, 1));  % momento angolare specifico: una volta plottato
% permette di capire il verso di rotazione
% dell'orbita

%momento della quantità di moto per disegno:

h=h/norm(h);  %normalizzo h in modo che la dimensione non influenzi la scala della rappresentazione 
h=3*R_rif*h;     % "allungo" h in modo che non venga coperto dalla superficie terrestre

hx=[0, h(1)];
hy=[0, h(2)];
hz=[0, h(3)];

%perigeo e apogeo per disegno:

xp=rr_min(1);  %x perigeo
yp=rr_min(2);  %y perigeo
zp=rr_min(3);  %z perigeo
if e<1
    xa=rr_max(1);  %x apogeo
    ya=rr_max(2);  %y apogeo
    za=rr_max(3);  %z apogeo
end

%plotting: 

if strlength(color)>1
    rz(c)=NaN;
    C=linspace(0, 100, c);
    p = patch(rx, ry, rz, C, 'FaceColor','none','EdgeColor', 'interp', 'linewidth', 1, DisplayName=OrbitName); %use patch   %'FaceColor','none','EdgeColor','interp'
    CB=colorbar();                                                      %'Direction', 'reverse'
    CB.Label.String = 'ToF% in the orbit';

    oldcmap = colormap(color);                                            %autumn(100)
    colormap( flipud(oldcmap) );

                                                                        %colormap(p, hot);    %autumn(100)
    hold on
    plot3(hx, hy, hz, 'k.-', HandleVisibility='off'); %Angular momentum DisplayName='Specific Angular Momentum'
else
    legend();
    width = 1;
    plot3(rx, ry, rz, '-', Color=color, DisplayName=OrbitName, LineWidth = width); %trajectory
    hold on
    plot3(hx, hy, hz, '.-', Color=color, HandleVisibility='off'); %momento della quantità di moto DisplayName='Starting Specific Angular Momentum'
end

quiver3(0,0,0,2*R_rif,0,0,'k','linewidth', 1, 'handlevisibility', 'off');   %versore I
quiver3(0,0,0,0,2*R_rif,0,'k', 'linewidth', 1, 'handlevisibility', 'off');  %versore J
quiver3(0,0,0,0,0,2*R_rif,'k', 'linewidth', 1 , 'handlevisibility', 'off'); %versore K

%aggiungo Marker in punto iniziale e finale

if th0==0 && thf==2*pi && nargin<15

else
    plot3(rr0(1),rr0(2), rr0(3), '.', 'MarkerSize', 16, MarkerEdgeColor='g', HandleVisibility='off');
    plot3(rrf(1),rrf(2), rrf(3), '.', 'MarkerSize', 16, MarkerEdgeColor='r', HandleVisibility='off');
end

%eventuale aggiunta di pericentro e apocentro

if nargin<9
        optn='null';
end

if optn=='n'
    plot3(xp, yp, zp, 'k.', HandleVisibility='off', LineWidth=18)
    if e<1          %if elliptic, plot the apocentre
        plot3(xa, ya, za, 'k.',  HandleVisibility='off', LineWidth=18);  
    elseif e>1      %if hyperbolic, plot the asymptots
        rp=sqrt(xp^2 + yp^2 + zp^2);
        k=(rp-a)/rp;
        Asy=zeros(3);
        Asy(2, :)=[xp*k, yp*k, zp*k];
        if th0*thf<0                 %show both the asymptot
            dh=asin(1/e); %half of turn angle
            L3=(ROT(i, Om, om)*rrmax*[cos(pi/2+dh), sin(pi/2+dh), 0]')';
            Asy(3, :)=Asy(2, :)+L3;

            L1=(ROT(i, Om, om)*rrmax*[cos(pi/2+dh), -sin(pi/2+dh), 0]')';
            Asy(1, :)=Asy(2, :)+L1;
        elseif th0*thf>=0 && max(th0, thf)>0    %show only the outgoing asymptot
            dh=asin(1/e); %half of turn angle
            L3=(ROT(i, Om, om)*rrmax*[cos(pi/2+dh), sin(pi/2+dh), 0]')';
            Asy(3, :)=Asy(2, :)+L3;

            Asy(1, :)=Asy(2, :);
        elseif th0*thf>=0 && min(th0, thf)<=0    %show only the ingoing asymptot
            dh=asin(1/e); %half of turn angle
            L1=(ROT(i, Om, om)*rrmax*[cos(pi/2+dh), -sin(pi/2+dh), 0]')';
            Asy(1, :)=Asy(2, :)+L1;

            Asy(3, :)=Asy(2, :);
        end

        if strlength(color)>1
            plot3(Asy(:, 1), Asy(:, 2), Asy(:, 3), 'k-.', HandleVisibility='off');
        else
            plot3(Asy(:, 1), Asy(:, 2), Asy(:, 3), '-.', HandleVisibility='off', Color=color);
        end
end

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title(Title);
axis equal;
grid on;

end


