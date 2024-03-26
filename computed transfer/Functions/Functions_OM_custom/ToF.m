function [T] = ToF(a, e, th0, thf, mu, p)
%ToF computes the time of flight given the semi-major axis, the
%eccentricity, the starting and final true anomaly, the gravitational
%constant and (necessary for parabolic orbit) the semi-parameter
%
%INPUTs:
% a   : semi-major axis            [1x1] [km]
% e   : eccentricity               [1x1] [-]
% th0 : starting true anomaly      [1x1] [rad]
% thf : final true anomaly         [1x1] [rad]
% mu  : gravitational constant     [1x1] [km^3/s^2]
% p   : semi-parameter             [1x1] [km] 
%OUTPUTs:
% T   : Time of Flight (ToF)       [1x1] [s]
%
% if thf<th0, the function gives a negative value for the ToF
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

if e==1 && a==Inf   %if parabolic orbit

    np=2*sqrt(mu/p^3);

    if nargin<6

        error("Inserire il valore di p = 2*rP negli Input!")

    end

    if th0==-pi  %controllo su esistenza di theta

        error("L'orbita parabolica non passa per -pi!")

    else
    
        B1=tan(th0/2);

        T1= (B1+(B1^3)/3) / np;

    end

    if thf==pi   %controllo su esistenza di theta

        error("L'orbita parabolica non passa per pi!")

    else
    
        B2=tan(thf/2);

        T2= (B2+(B2^3)/3) / np;

    end

    T=T2-T1;

elseif e<1 && e>=0     %se orbita ellittica

    n=sqrt(mu/a^3);  %velocità angolare media

    Tp=2*pi/n;   % periodo orbitale nominale

    if th0==0
        th1=0;
    else
        th1=mod(abs(th0), 2*pi)*(th0/abs(th0)); %resto del numero di periodi completati da th0
    end
    K1=(th0-th1)/(2*pi);  %numero di periodi completati da th0

    if thf==0
        th2=0;
    else
        th2=mod(abs(thf), 2*pi)*(thf/abs(thf)); %resto del numero di periodi completati da thf
    end
    K2=(thf-th2)/(2*pi);  %numero di periodi completati da thf


    S1= sin(th1) * sqrt(1-e^2) / (1 + e*cos(th1)); %sin(E1)
    C1=         (e + cos(th1)) / (1 + e*cos(th1)); %cos(E1)
    if atan2(S1, C1)>0 && th1>0
        E1= atan2(S1, C1);
    elseif atan2(S1, C1)>0 && th1<0
        E1= atan2(S1, C1) - 2*pi;
    elseif th1>0 % && atan2<0
        E1= atan2(S1, C1) + 2*pi;
    else %th1<0 && atan2<0
        E1= atan2(S1, C1);
    end
    t1= (E1 - e*S1) / n; %tempo per andare da theta=0 a th1


    S2= sin(th2) * sqrt(1-e^2) / (1 + e*cos(th2)); %sin(E2)
    C2=         (e + cos(th2)) / (1 + e*cos(th2)); %cos(E2)
    if atan2(S2, C2)>0 && th2>0
        E2= atan2(S2, C2);
    elseif atan2(S2, C2)>0 && th2<0
        E2= atan2(S2, C2) - 2*pi;
    elseif th2>0 % && atan2<0
        E2= atan2(S2, C2) + 2*pi;
    else %th2<0 && atan2<0
        E2= atan2(S2, C2);
    end
    t2= (E2 - e*S2) / n; %tempo per andare da theta=0 a th2

    
    T1= K1*Tp +t1;   %tempo per andare da theta=0 a th0

    T2= K2*Tp +t2;   %tempo per andare da theta=0 a thf

    T= T2 - T1;  %tempo complessivo per andare da th0 a thf

elseif e>1        %se orbita iperbolica

    if a>=0  %controllo su input
        error("Il semiasse maggiore di un'orbita iperbolica deve essere negativo in segno!")
    end

    alpha=pi-acos(1/e);      %anomalia limite per orbita iperbolica: 2*alpha è l'angolo (maggiore) tra i due asintoti

    if th0<=-alpha   %controllo su esistenza di theta
        warning("L'orbita parabolica non passa per anomalie inferiori all'angolo limite!")
    elseif th0>=alpha   %controllo su esistenza di theta
        warning("L'orbita parabolica non passa per anomalie inferiori all'angolo limite!")
    else
    
        F1=2*atanh(tan(th0/2)*sqrt((e-1)/(e+1)));

        tan(th0/2)*sqrt((e-1)/(e+1));

        T1=sqrt(((-a)^3)/mu)*(e*sinh(F1)-F1);

    end

    if thf>=alpha   %controllo su esistenza di theta
        error("L'orbita parabolica non passa per anomalie superiori all'angolo limite!")
    elseif thf<=-alpha %controllo su esistenza di theta
        error("L'orbita parabolica non passa per anomalie superiori all'angolo limite!")
    else
    
        F2=2*atanh(tan(thf/2)*sqrt((e-1)/(e+1)));

        tan(thf/2)*sqrt((e-1)/(e+1));

        T2=sqrt(((-a)^3)/mu)*(e*sinh(F2)-F2);

    end

    T=T2-T1;

end
end

