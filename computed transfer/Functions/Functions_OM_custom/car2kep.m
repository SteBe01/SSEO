function [a, e, i, OM, om, f, p] = car2kep(rr, vv, mu)
%Given the position vector rr, the velocity vector vv and the gravitational
%costant mu, calculates the keplerian coordinates of the orbit
%
%INPUTS:
% rr : position vector         [3x1] [km]
% vv : velocity vector         [3x1] [km/s]
% mu : gravitational constant  [1x1] [km^3/s^2]
%
%OUTPUTS:
% a  : semi-major axis                                                                                     [1x1] [km]
% e  : eccentricity                                                                                        [1x1] [-]
% i  : inclination        ( 0=< i =<pi )                                                                   [1x1] [rad]
% OM : right ascension of the ascending node (RAAN) ( 0=<OM=<2*pi )                                        [1x1] [rad]
% om : pericenter anomaly ( 0=<om=<2*pi )                                                                  [1x1] [rad]
% f  : true anomaly       ( 0=<f=<2*pi )                                                                   [1x1] [rad]
% p  : semi-parameter (distance from the focus when f=pi/2: needed when handling a parabolic orbit)        [1x1] [km]
%
% this function handles values expressed in the IS (km, s, etc.)
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.
%
%debugged!
%
%previously known as "rv2parorb_p"
%M=theta=f

if size(rr)==[1, 3]  %the function works with column vectors only!
    rr=rr';
end

if size(vv)==[1, 3]  %the function works with column vectors only!
    vv=vv';
end

r=norm(rr); %norma posizione
v=norm(vv); %norma velocità
eps=((v^2)/2)-(mu/r); %energia meccanica specifica
vr=dot(rr, vv)/r; %velocità radiale


% controllo su forma orbita usando eps

a=-mu*0.5/eps; %semiasse maggiore

ee=((v^2-mu/r)*rr-(r*vr)*vv)/mu; %vettore eccentricità
e=norm(ee); %eccentricità
h=cross(rr, vv); %momento angolare
i=acos(h(3)/norm(h));
k=[0 0 1]'; %versore z
N=cross(k, h); %linea dei nodi

if N(2)>0
    OM=acos(N(1)/norm(N));
elseif N(2)<0
    OM=2*pi-acos(N(1)/norm(N));
elseif N(2)==0
    OM=0;
end


if ee(3)>0
    om=acos(dot(N, ee)/(norm(N)*e));
elseif ee(3)<0
    om=2*pi-acos(dot(N, ee)/(norm(N)*e));
elseif ee(3)==0
    om=0;
end

if eps<0
    if vr>0
        f=acos(dot(ee, rr)/(e*r));
    elseif vr<0
        f=2*pi-acos(dot(ee, rr)/(e*r));
    elseif vr==0
        f=0;
    end
elseif eps>=0
    if vr>0
        f=acos(dot(ee, rr)/(e*r)) ;
    elseif vr<0
        f=-acos(dot(ee, rr)/(e*r)) ;
    elseif vr==0
        f=0;
    end
end

p=r*(1+e*cos(f));  %valid for every kind of orbit

if e>=1 && a>0                                  %controllo eccentricità

    disp("l'orbita è parabolica");

elseif e>1                             %controllo eccentricità

    disp("l'orbita è iperbolica");

elseif e<1 && e>0                      %controllo eccentricità

    %disp("l'orbita è ellittica");
    
elseif e==0                            %controllo eccentricità

    %disp("l'orbita è circolare");    

end

