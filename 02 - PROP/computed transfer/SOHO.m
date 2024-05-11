%%

clc, clear
close all

restoredefaultpath
addpath(genpath('Functions'))

mu_Sun = astroConstants(4); %[km^3/s^2] Sun gravity constant
mu_Earth = astroConstants(13); %[km^3/s^2] Earth gravity constant
au = astroConstants(2); %[km] Astronomic Unit
r_Earth=astroConstants(23);%[km] mean radius of the Earth
lCCAS=deg2rad(28.4555555); %[rad] Cape Canaveral Latitude;
w_E=15.04*pi/(180*3600);%[rad/s] Earth spin rate

r_INJECTION=1380360; %[km] distance from the Earth at which 
Ax = 206448; %[km] Semi-Diameter of "orbit" around L1
Ay = 666672; %[km] Semi-Diameter of "orbit" around L1
Az = 120000; %[km] Semi-Diameter of "orbit" around L1


% Coordinates of Earth at departure of the mission:
[kEa1,   ~] = uplanet(date2mjd2000([1995, 12, 2, 8, 8, 1]), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa1, vEa1] = kep2car(kEa1(1), kEa1(2), kEa1(3), kEa1(4), kEa1(5), kEa1(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
rpEa = kEa1(1)*(1-kEa1(2)); %pericentre radius of Earth
nEa  = sqrt(mu_Sun/kEa1(1)^3);  %medium angular velocity of the Earth
r_SUN=norm(rEa1); %[km] distance of Earth wrt Sun in 2/12/1995 
rSOI1=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth


% Coordinates of Earth at injection in L1 orbit:
[kEa2,   ~] = uplanet(date2mjd2000([1996, 2, 14, 0, 0, 0]), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
r_SUN=norm(rEa2); %[km] distance of Earth wrt Sun in 14/2/1996 
rSOI2=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth

%---> SOHO in the ending orbit (L1) is out of Earth SOI

%DeltaV to escape Earth from Cape Canaveral
ve0=sqrt(2*mu_Earth/r_Earth);%[km/s] overall deltaV necessary to escape Earth Gravity
v0=cos(lCCAS)*r_Earth*w_E;%[km/s] Cape Canaveral tangential velocity due to Earth spin rate
dv0=ve0-v0;%[km/s] deltaV expense to escape Earth Gravity

theta_SOI = acos((r_Earth/max(rSOI1, rSOI2))-1);
ToF_SOI_parabolic=ToF(Inf, 1, pi/2, theta_SOI, mu_Earth, r_Earth);
ToF_SOI_parabolic_days=ToF_SOI_parabolic/(24*3600);

theta_INJECTION = acos((r_Earth/r_INJECTION)-1);
ToF_INJECTION=ToF(Inf, 1, pi/2, theta_INJECTION, mu_Earth, r_Earth);
ToF_INJECTION_days=ToF_INJECTION/(24*3600);

% Position of L1 at departure date:
f = @(b) 1+(mu_Earth/mu_Sun)/((1-b)*b^2)-1/(1-b)^3 ;

b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)

rLag1=rEa1*(1-b);%[km] L1 position wrt to Sun on 2/12/1995
vLag1=vEa1*(1-b);%[km/s] L1 velocity wrt to Sun on 2/12/1995

% Position of L1 at injection date:
rLag2=rEa2*(1-b);%[km] L1 position wrt to Sun on 14/2/1996
vLag2=vEa2*(1-b);%[km/s] L1 velocity wrt to Sun on 14/2/1996

distance_ratio1=norm(b*rEa1)/rSOI1; %ratio between distance of Earth from L1 and Earth SOI radius
distance_ratio2=norm(b*rEa2)/rSOI2; %ratio between distance of Earth from L1 and Earth SOI radius

[t, r, v] = extract_horizons('soho_ins_fine.txt');

n=size(t, 2);

for i=1:n-1
    dv_INJ(i)=norm(v(:, i+1)-v(:, i));
end
dv_INJ(n)=0;

plot(t, dv_INJ)

% %Use Lambert to compute the extimated espense
% %parameters of Lambert function:
% orbitType = 0;    %prograde or retrograde (clockwise=0, counterclockwise=1)
% Nrev = 0 ;        %number of revolutions (Nrev = 0 ZERO-REVOLUTION transfer is calculated, Nrev > 0 two transfers are possible)
% Ncase = 0 ;       %case of a (in case of Nrev>0:    0: small-a option, 1: large-a option)
% optionsLMR = 1;   %select kind of errors shown by Lambert
% ToF=74 * 24 * 3600; %[s] time of Flight
% 
% %First estimation (departure at center of Earth, Arrival in L1, tof=74
% %days)
% [A1,P1,E1,ERROR1,VI1,VF1,TPAR1,THETA1] = lambertMR(rEa1,rLag2,ToF,mu_Sun,orbitType,Nrev,Ncase,optionsLMR);
% dvI1=norm(VI1-vEa1);
% dvF1=norm(VF1-vLag2);
% dv1=dvI1+dvF1
% 
% PlotOrbitWtime(rEa1, vEa1, ToF, mu_Sun, 'g', 'SOHO', 'Earth');
% 
% PlotOrbitWtime(rEa1, VI1, ToF, mu_Sun, 'm', 'SOHO', 'Transfer');
% %PlotOrbitWtime(rLag1, vLag1, ToF, mu_Sun, 'r', 'SOHO', 'L1');

% r1  = NaN*ones(3,l1) ;
% dV1 = NaN*ones(l1,m1) ;  
% Ta1 = NaN*ones(l1,m1) ; 
% Tp1 = NaN*ones(l1,m1) ;
% Te1 = NaN*ones(l1,m1) ;
% Ter1 = NaN*ones(l1,m1) ;
% Tv1 = NaN*ones(l1,m1, 3) ;
% Tv2 = NaN*ones(l1,m1, 3) ;
% Ttpa1 = NaN*ones(l1,m1, 3) ;
% Tth1 = NaN*ones(l1,m1, 3) ;
% tof1 = NaN*ones(l1,m1, 3) ;
% v_Ea1= NaN*ones(l1,m1, 3) ;
% 
% for i = 1:l1
%     [kep1, ~] = uplanet(tspan1(i), 3) ;   %get the keplearian coordinates of Earth  
%     [r1(:, i),v1] = kep2car(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), mu_Sun) ; %get the position and velocity of Earth
% 
%     for j = 1:m1
%         if tspan2(j)<=tspan1(i)
%             continue
%         else
%             [kep2, ~] = uplanet(tspan2(j), ID2) ; %get the keplerian coordinates of Saturn
%             [r2,v2] = kep2car(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), mu_Sun) ; %get the position and velocity of Earth
% 
%             tof1(i, j) = (tspan2(j) - tspan1(i))*86400 ;   %Time of Flight in seconds
%             [Ta1(i, j),Tp1(i, j),Te1(i, j),Ter1(i, j),Tv1(i, j, :),Tv2(i, j, :),Ttpa1(i, j),Tth1(i, j)] = lambertMR(r1(:, i),r2,tof1(i, j),mu_sun,orbitType,Nrev,Ncase,optionsLMR) ;
% 
%             v_Ea1(i, j, :)=[Tv2(i, j, 1)-v2(1), Tv2(i, j, 2)-v2(2), Tv2(i, j, 3)-v2(3)]';
%             n1=sqrt((Tv1(i, j, 1) - v1(1))^2 + (Tv1(i, j, 2) - v1(2))^2 + (Tv1(i, j, 3) - v1(3))^2 );
%             n2=sqrt((Tv2(i, j, 1) - v2(1))^2 + (Tv2(i, j, 2) - v2(2))^2 + (Tv2(i, j, 3) - v2(3))^2 );
% 
%             dV1(i,j) = n1; 
%         end
%     end
% end


% h1 = 0.1; %time discretization [days]
% 
% % Defining possible departure dates
% fdd1 = date2mjd2000([1995, 12, 1, 0, 0, 0]); %first possible departure day
% ldd1 = date2mjd2000([1995, 12, 3, 0, 0, 0]); %last possible departure day
% 
% tspan1 = fdd1 : h1 : ldd1 ;
% 
% % Defining possible arrival dates: 
% fad1 = date2mjd2000([1996, 2, 13, 0, 0, 0]);  %first possible arrival day
% lad1 = date2mjd2000([1996, 2, 15, 0, 0, 0]);  % last possible arrival day
% 
% tspan2 = fad1 : h1 : lad1; 
% 
% l1=length(tspan1);
% m1=length(tspan2);

% %orbital parameters around L1 point
% 
% a = 148.110e6; %[km] semimajor axis of the ending orbit in Sun frame
% 
% T_L1 = 178 * 84600; %[s] approximate orbital period about L1
% 
% sd1 = Ay; %[km] first semidiameter of the orbit around L1
% sd2 = sqrt(Ax^2 + Az^2); %[km] second semidiameter of the orbit around L1
% 
% % the insertion takes place in the 10th week from the starting of the
% % mission
% 
% ra = sqrt( (a+Ax)^2 + Az^2 ); %[km] apocenter radius of the ending orbit in Sun frame
% rp = sqrt( (a-Ax)^2 + Az^2 ); %[km] pericenter radius of the ending orbit in Sun frame
% 
% e = (ra-rp)/(ra+rp); %[-] eccentricity of the ending orbit in Sun frame
% 
% i = atan2(Az, a+Ax); %[rad] inclination of ending orbit in Sun frame
% i_deg = rad2deg( i ); %[deg] inclination of ending orbit in Sun frame
% 
% kL1=kEa1;
% kL1(1)=a;
% kL1(2)=e;
% kL1(3)=i;
% %kL1(4)=...  (OM)
% %kL1(5)=...  (om)
% %kL1(6)=...  (theta)