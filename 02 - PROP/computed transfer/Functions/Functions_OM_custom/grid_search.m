function [dv, dv1, dv2, dv3, d, rp, tdep, tGA, tarr] = grid_search(tspan,dt,h,ID1,ID2,ID3,orbitType,Nrev,Ncase,optionsLMR)
%GRID_SEARCH computes the many different values of dv, dv1, dv2, dv3, d,
%rp, tdep, tGA, tarr in a radius of dt, with a discretization h.
%
% INPUTS:
% tspan  : vector of dates for the mission (dep, ga, arr) in mjd2000            [3x1]   [days]
% dt     : radius of the neighbourhood of research for each date (dep, ga, arr) [1x1]   [days]
% h      : vector of the discretizations of the neighbourhoods of each date     [3x1]   [days]
% ID1    : ID of the departure celestial body                                   [1x1]   [-]
% ID2    : ID of the flyby celestial body                                       [1x1]   [-]
% ID3    : ID of the arrival celestial body                                     [1x1]   [-]
%
%	orbitType[1]    Logical variable defining whether transfer is
%                       0: direct transfer from R1 to R2 (counterclockwise)
%                       1: retrograde transfer from R1 to R2 (clockwise)
%	Nrev[1]         Number of revolutions.
%                   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%                   if Nrev > 0 two transfers are possible. Ncase should be
%                          defined to select one of the two.
%	Ncase[1]        Logical variable defining the small-a or large-a option
%                   in case of Nrev>0:
%                       0: small-a option
%                       1: large-a option
%	optionsLMR[1]	lambertMR options:
%                    optionsLMR(1) = display options:
%                                    0: no display
%                                    1: warnings are displayed only when
%                                       the algorithm does not converge
%                                    2: full warnings displayed
%
% OUTPUTS:
% dv  : 3-D tensor storing the total dv for each mission                  [l x m x n]   [km/s]
% dv1 : matrix storing the departure expense for each mission             [l x m ]      [km/s]
% dv2 : 3-D tensor storing the gravity assist expense for each mission    [l x m x n]   [km/s]
% dv3 : 3-D tensor storing the arrival expense for each mission           [l x m x n]   [km/s]
% d   : 3-D tensor storing the turn angle of flyby of each mission        [l x m x n]   [rad]
% rp  : 3-D tensor storing the pericentre of flyby of each mission        [l x m x n]   [km]
% tdep: vector of tested dates for departure in mjd2000                   [l x 1]       [days]
% tGA : vector of tested dates for gravity assist in mjd2000              [m x 1]       [days]
% tarr: vector of tested dates for arrival in mjd2000                     [n x 1]       [days]
%
% if tGA(j) is less than tdep(i), or if tarr(k) is less than tGA(j), the
% function returns dv(i, j, k)=NaN, and so o for dv1, dv2, dv3, d, rp
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

if size(tspan)~=[3, 1]
    tspan=tspan';
end

h1=h(1);
h2=h(2);
h3=h(3);

l=1+floor(2*dt/h1);
m=1+floor(2*dt/h2);
n=1+floor(2*dt/h3 );

% Defining ranges of possible dates for departure, flyby and arrival
tlim1=date2mjd2000([2028, 1, 1, 0, 0, 0]);
tlim2=date2mjd2000([2058, 1, 1, 0, 0, 0]);

start_dep = tspan(1) - dt ; 
final_dep = tspan(1) + dt ; 

start_GA  = tspan(2) - dt ; 
final_GA  = tspan(2) + dt ; 

start_arr = tspan(3) - dt ; 
final_arr = tspan(3) + dt ;


tdep = linspace( max(start_dep, tlim1), min(final_dep, tlim2), l) ;
tGA  = linspace( max(start_GA,  tlim1), min(final_GA,  tlim2), m) ; 
tarr = linspace( max(start_arr, tlim1), min(final_arr, tlim2), n) ;

% constants:
mu_sun = astroConstants(4) ;  

dv1=ones(l, m)*NaN ;
dv2=ones(l, m, n)*NaN ;
rp =ones(l, m, n)*NaN ;
dv3=ones(l, m, n)*NaN ;
dv =ones(l, m, n)*NaN ;
d  =ones(l, m, n)*NaN ;

for i = 1:l
    if ID1<11
        [kSa, ~] = uplanet(tdep(i), ID1) ;   %get the keplearian coordinates of Saturn 
    else
        [kSa, ~] = ephNEO(tdep(i), ID1) ;   %get the keplearian coordinates of Saturn
    end
    [rSa,vSa] = kep2car(kSa(1), kSa(2), kSa(3), kSa(4), kSa(5), kSa(6), mu_sun) ; %get the position and velocity of Saturn

    for j = 1:m
        if tGA(j)<=tdep(i) + 2044/2 %half of the minimum time for a bitangent transfer
            continue
        elseif tGA(j)>=tdep(i) + 3*2044 %3 times of the minimum time for a bitangent transfer
            continue
        else
            if ID2<11
                [kEa, ~] = uplanet(tGA(j), ID2) ; %get the keplerian coordinates of the Earth
            else
                [kEa, ~] = ephNEO(tGA(j), ID2) ; %get the keplerian coordinates of the Earth
            end
            [rEa,vEa] = kep2car(kEa(1), kEa(2), kEa(3), kEa(4), kEa(5), kEa(6), mu_sun) ; %get the position and velocity of Earth
        
            tof1 = (tGA(j) - tdep(i))*86400 ;   %Time of Flight 1
            [~,~,~,~,v1,v2,~,~] = lambertMR(rSa,rEa,tof1,mu_sun,orbitType,Nrev,Ncase,optionsLMR) ;
            dv1(i, j)= norm(v1' - vSa);
        end

        for k = 1:n
            if tarr(k)<=tGA(j)
                continue
            else
                if ID3<11
                    [k65, ~] = uplanet(tarr(k) , ID3) ; %get the keplerian coordinates of asteroid 65
                else
                    [k65, ~] = ephNEO(tarr(k) , ID3) ; %get the keplerian coordinates of asteroid 65
                end
                [r65,v65] = kep2car(k65(1), k65(2), k65(3), k65(4), k65(5), k65(6), mu_sun) ; %get the position and velocity of asteroid 65

                tof2 = (tarr(k) - tGA(j))*86400 ;   %Time of Flight 2
                [~,~,~,~,v3,v4,~,~] = lambertMR(rEa,r65,tof2,mu_sun,orbitType,Nrev,Ncase,optionsLMR) ;

                vI=v2'-vEa;
                vF=v3'-vEa;

                d(i, j, k)=real(acos(dot(vI, vF)/(norm(vI)*norm(vF)))) ;   %delta
            
                [dv2(i, j, k), rp(i, j, k), ~,~]=fly_by(vI, vF, ID2);   %dv for flyby

                dv3(i, j, k)=norm(v4'-v65);
         
                dv(i,j, k) =dv1(i, j)+dv2(i, j, k)+dv3(i, j, k);
            end
        end
    end
end




end

