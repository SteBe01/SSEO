function [th_0, th_end] = PaleBluePlot(Planet_ID, T0, Tend, color, Title, OrbitName, ~)
%GOTOIN Summary of this function goes here
%THIS FUNCTION PLOTS THE PLANET ORBIT (AROUND THE SUN!!) GIVEN THE STARTING
%DATE AND THE ARRIVAL DATE, AND GIVES AS OUTPUTS THE STARTING AND FINAL
%TRUE ANOMALIES
% INPUTS:
%	Planet_ID    Integer number identifying the celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   
%   T0:    starting date in form of vector [year, month, day, hour, minute, second]
%
%   Tend:   arrival date in form of vector [year, month, day, hour, minute, second]
%
%   ~: options: if any input, the function will also plot the pericentre and
%            apocentre, otherwise the function will show only the orbit
%
% OUTPUTS:
% th_0    : initial true anomaly [rad]   [1x1]
% th_end  : finale true anomaly  [rad]   [1x1]
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

%
%The name of this function was inspired by Carl Sagan's speech "Pale Blue Dot"
%

mjd2000_0= date2mjd2000(T0); %starting date in days

mjd2000_end= date2mjd2000(Tend); %arrival date in days

T=(mjd2000_end-mjd2000_0)*24*3600; %ToF in seconds

%get the keplerian coordinates of the planet:

[kep_0,   ~]=uplanet(mjd2000_0, Planet_ID);

[kep_end, ~]=uplanet(mjd2000_end, Planet_ID);

mu_sun=astroConstants(4);

%extract the true anomalies:

th_0  =kep_0(6); %starting true anomaly

th_end=kep_end(6); %arrival true anomaly

%plot the orbit:
a =kep_0(1);  %semimajor axis

e =kep_0(2);  %eccentricity

i =kep_0(3);  %inclination

Om=kep_0(4);  %RAAN (right ascension of the ascending node)

om=kep_0(5);  %pericentre anomaly

%get the initial position and velocity:

[r0, v0]=kep2car(a, e, i, Om, om, th_0, mu_sun);

% if nargin>6 %show apocentre and pericentre  (~=options='n')
%     %plotOrbit(a, e, i, Om, om, th_0, th_end, mu_sun, [], 'n', color, Title, OrbitName);
% else %omit apocentre and pericentre
%     %plotOrbit(a, e, i, Om, om, th_0, th_end, mu_sun, [], 'no', color, Title, OrbitName);
% end

p=floor(ToF(a, e, 0, 2*pi, mu_sun)/(365.25*86400)); %approximate period of the celestial body

if Tend(1)==T0(1)+p
    plotOrbit(a, e, i, Om, om, 0, 2*pi, mu_sun, [], 'no', color, Title, OrbitName);
else
    PlotOrbitWtime(r0, v0, T, mu_sun, color, Title, OrbitName);
end


end

