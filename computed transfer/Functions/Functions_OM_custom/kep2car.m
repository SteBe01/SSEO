function [rr, vv] = kep2car(a, e, i, OM, om, f, mu, p)
%Given the keplerian coordinates computes the cartesian coordinates
%(position and velocity vectors)
%
%INPUTS:
% a  : semi-major axis                                                                               [1x1] [km]
% e  : eccentricity                                                                                  [1x1] [-]
% i  : inclination ( 0=< i =<pi )                                                                    [1x1] [rad] 
% OM : right ascension of the ascending node (RAAN) ( 0=<OM=<2*pi )                                  [1x1] [rad] 
% om : pericenter anomaly ( 0=<om=<2*pi )                                                            [1x1] [rad] 
% f  : true anomaly ( 0=<f=<2*pi )                                                                   [1x1] [rad]   
% mu : gravitational constant                                                                        [1x1] [km^3/s^2]
% p  : semi-parameter (distance from the focus when f=pi/2: needed when handling a parabolic orbit)  [1x1] [km]
%
%
%OUTPUTS:
% rr : position vector [3x1] [km]
% vv : velocity vector [3x1] [km/s]
%
% this function handles values expressed in the IS (km, s, etc.)
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.
%
%debugged!
%
%previously known as "parorb2rv_p"
%
%M=theta=f

if i<0
    warning("i cannot be less than zero!");
    i=abs(mod(i, pi));
elseif i>pi
    warning("i cannot be more than pi!");
    i=mod(i, pi);
elseif e<0
    warning("e cannot be less than zero!");
elseif OM<0
    OM=mod(OM, 2*pi)+2*pi;
elseif OM>2*pi
    OM=mod(OM, 2*pi);
elseif om<0
    om=mod(om, 2*pi)+2*pi;
elseif om>2*pi
    om=mod(om, 2*pi);
end

if e~=1

    if e>1                      %controllo eccentricità

        if a>=0

            error("Il semiasse maggiore di un'orbita iperbolica deve essere negativo in segno!");

        end

        %p=a*(1-e)*(1+e);     %semilato retto
        
        alpha=pi-acos(1/e);    %anomalia limite per orbita iperbolica


        if f<=-alpha     %controllo su anomalia massima se orbita iperbolica
             
            error("f non appartiene all'orbita iperbolica!");

        end

        if f>=alpha     %controllo su anomalia massima se orbita iperbolica
             
            error("f non appartiene all'orbita iperbolica!");

        end

    end

p=a*(1-e)*(1+e); %semilato retto

r=p/(1+e*cos(f));                  %modulo raggio

v_theta=((mu/p)^0.5)*(1+e*cos(f));  %modulo velocità tangenziale
v_r=((mu/p)^0.5)*(e*sin(f));        %modulo velocità radiale

rr=[r*cos(f); r*sin(f); 0]; %vettore posizione su piano xy (colonna)
v_rad=rr*(v_r/r);                   %vettore velocità radiale in piano xy
%if i<=pi/2;
    k=[0 0 1]'; %versore parallelo a z
v_tan=v_theta*(cross(k, rr)/r);   %vettore velocità tangenziale in piano xy.
                                  %ho dei dubbi sul calcolo del verso perché questo dipende da
                                  %i, ma questo entra in gioco quando al
                                  %termine del processo ruoto tutto il
                                  %sistema di riferimento costruito.
                                  
vv=[v_rad+v_tan];                   %vettore velocità in piano xy (colonna)

elseif e==1 && a>1e20

    if nargin<8

        error("Inserire il valore di p = 2*rP negli Input!");

    end

    alpha=pi;   %anomalia limite per orbite paraboliche
   
    if f<=-alpha

        error("f non appartiene all'orbita parabolica!");

    end

    if f>=alpha

        error("f non appartiene all'orbita parabolica!");

    end

    r=p/(1+cos(f));                  %modulo raggio

    v_theta= sqrt(mu/p) * (1+cos(f));  %modulo velocità tangenziale

    v_r=     sqrt(mu/p) * (sin(f));    %modulo velocità radiale

    rr=[r*cos(f); r*sin(f); 0]; %vettore posizione su piano xy (colonna)

    v_rad=v_r*(rr/r);                   %vettore velocità radiale in piano xy

    k=[0 0 1]';                       %versore parallelo a z

    v_tan=v_theta*(cross(k, rr)/r);   %vettore velocità tangenziale in piano xy.
                                      %ho dei dubbi sul calcolo del verso perché questo dipende da
                                      %i, ma questo entra in gioco quando al
                                      %termine del processo ruoto tutto il
                                      %sistema di riferimento costruito.
                                  

    vv=[v_rad+v_tan];                 %vettore velocità in piano xy (colonna)

end

R=ROT(i, OM, om);   %prodotto del trasposto del prodotto delle matrici di rotazione (matrice 3x3)

rr=R*rr; %vettore posizione ruotato in piano orbitale (colonna)
vv=R*vv; %vettore velocità ruotato in piano orbitale  (colonna)

end