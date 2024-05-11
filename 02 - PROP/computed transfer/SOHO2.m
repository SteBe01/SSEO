%%
clc;
clear;
close all;
addpath(genpath('Functions'));

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
k=1;
for i=1:n-1
    dv_INJ(i)=norm(v(:, i+1)-v(:, i));
    if i>1
        if abs(dv_INJ(i))>= abs(1.01*dv_INJ(i-1))
            pos(k)=i;
            k=k+1;
        elseif abs(dv_INJ(i))<= abs(0.99*dv_INJ(i-1))
            pos(k)=i;
            k=k+1;
        end
    end
end

for i=pos(1):pos(end)
    dv_IST(i)=dv_INJ(i)-dv_INJ(pos(1)-1);
end
dv_INJEC=sum(dv_IST);

dv_INJ(n)=0;
figure(1)
plot(t, dv_INJ);

INJ_manoeuvre_duration=t(pos(end))-t(pos(1));

date_starting_INJ=jd2date(t(pos(1)));
date_ending_INJ=jd2date(t(pos(end)));

dv_INJECTION=norm( (v(:, pos(1))-v(:, pos(1)-1)) - (v(:, pos(end))-v(:, pos(end)+1)));

[t_r, r_r, v_r] = extract_horizons('horizons_results-1.txt');

n_r=size(t_r, 2);
k=1;
for i=1:n_r-1
    dv_INJ_r(i)=norm(v_r(:, i+1)-v_r(:, i));
    if i>1
        if abs(dv_INJ_r(i))>= abs(1.01*dv_INJ_r(i-1))
            pos_r(k)=i;
            k=k+1;
        elseif abs(dv_INJ_r(i))<= abs(0.99*dv_INJ_r(i-1))
            pos_r(k)=i;
            k=k+1;
        end
    end
end

dv_INJ_r(n_r)=0;
figure(3)
plot(t_r, dv_INJ_r);

for i=2:135
    dv_r_1(i)=dv_INJ_r(i)-dv_INJ_r(1);
end
date_r_1_starting=jd2date(t_r(2));
date_r_1_ending=jd2date(t_r(135));
dv__r_1=sum(dv_r_1);

[t_l, r_l, v_l] = extract_horizons('horizons_results.txt');

n_l=size(t_l, 2);
k=1;
for i=1:n_l-1
    dv_INJ_l(i)=norm(v_l(:, i+1)-v_l(:, i));
    if i>1
        if abs(dv_INJ_l(i))>= abs(1.01*dv_INJ_l(i-1))
            pos_l(k)=i;
            k=k+1;
        elseif abs(dv_INJ_l(i))<= abs(0.99*dv_INJ_l(i-1))
            pos_l(k)=i;
            k=k+1;
        end
    end
end

for i=193:208
    dv_l_1(i)=dv_INJ_l(i)-dv_INJ_l(192);
end
date_l_1_starting=jd2date(t_l(193));
date_l_1_ending=jd2date(t_l(208));
for i=336
    dv_l_2(i)=dv_INJ_l(i)-dv_INJ_l(335);
end
date_l_2_starting=jd2date(t_l(336));
date_l_2_ending=jd2date(t_l(336));
for i=805:806
    dv_l_3(i)=dv_INJ_l(i)-dv_INJ_l(804);
end
date_l_3_starting=jd2date(t_l(805));
date_l_3_ending=jd2date(t_l(806));
for i=2146:2149
    dv_l_4(i)=dv_INJ_l(i)-dv_INJ_l(2145);
end
date_l_4_starting=jd2date(t_l(2146));
date_l_4_ending=jd2date(t_l(2149));
for i=3840:3842
    dv_l_5(i)=dv_INJ_l(i)-dv_INJ_l(3839);
end
date_l_5_starting=jd2date(t_l(3840));
date_l_5_ending=jd2date(t_l(3842));

dv__l_1=sum(dv_l_1);
dv__l_2=sum(dv_l_2);
dv__l_3=sum(dv_l_3);
dv__l_4=sum(dv_l_4);
dv__l_5=sum(dv_l_5);

dv_INJ_l(n_l)=0;
figure(2)
plot(t_l, dv_INJ_l);

dv_tot   =(dv0+dv__l_1+dv__l_2+dv__l_3+dv__l_4+dv__l_5)
dv_tot_sm=(dv0+dv__l_1+dv__l_2+dv__l_3+dv__l_4+dv__l_5)*1.1

%% Stima del dv totale usando integrali
clc; 
clear;
close all;
T_TC=478309; %[kN] spinta di ogni Thiokol Castor
T_AM=2093.7; %[kN] spinta di Atlas AM-5AS
T_AS=386300; %[kN] spinta di Atlas II-AS
T_CA=185012; %[kN] spinta di Ceantaur II-A

T0=4*T_TC+T_AM+T_AS;%[kN] spinta STADIO 0
T1=T_AM+T_AS;%[kN] spinta STADIO 1
T2=T_CA;%[kN] spinta STADIO 2

t_bTC=56; %[s] burning time Thiokol Castor
t_bAM=283; %[s] burning time Atlas AM-5AS
t_bAS=283; %[s] burning time Atlas II-AS
t_bCA=392; %[s] burning time Ceantaur II-A

t0=t_bTC; %[s] burning time stage 0
t1=t_bAM-t_bTC; %[s] burning time stage 1
t2=t_bCA; %[s] burning time stage 2

M0=235477; %[kg] mass at launch
M1=163771.53; %[kg] mass at ending stage 0
M2=157655.53; %[kg] mass at detach stage 0
M3=32605; %[kg] mass at ending stage 1
M4=20923; %[kg] mass at detach stage 1
M5=4143; %[kg] mass at ending stage 2
M6=1850; %[kg] mass at detach stage 2
M7=1599; %[kg] mass at ending SOHO tank

Q0=(M0-M1)/t0;
Q1=(M2-M3)/t1;
Q2=(M4-M5)/t2;
%Q3=(M6-M7)/


f0 = @(t) (T0)./(M0-t*Q0);
DV0=integral(f0, 0, t0);
f1 = @(t) (T1)./(M2-t*Q1);
DV1=integral(f1, 0, t1);
f2 = @(t) (T2)./(M4-t*Q2);
DV2=integral(f2, 0, t2);

t=[0:t0];
figure(1)
plot(t, f0(t), 'g')
t=[0:t1];
hold on
plot(t, f1(t), 'y')
t=[0:t2];
hold on
plot(t, f2(t), 'r')

DV3=220*9.81*log(M6/M7);

dvlauncher=DV0+DV1+DV2+DV3;

DV_tot=dvlauncher/1000; %é minore del valore per compiere la missione... male male


%% SOHO MISSION SIMULATION WITH THREE BODY PROBLEM SOLVING


clc; clear;
close all;
addpath(genpath('Default_Functions\'));
addpath(genpath('timeConversion\'));

G=6.6743e-20; %[km^3/(kg*s)]

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

Ax = 200000;
Ay = 650000;
Az = 200000;


% Coordinates of Earth at departure of the mission:
[kEa1,   ~] = uplanet(date2mjd2000([1995, 12, 2, 8, 8, 1]), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa1, vEa1] = kep2car(kEa1(1), kEa1(2), kEa1(3), kEa1(4), kEa1(5), kEa1(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
rpEa = kEa1(1)*(1-kEa1(2)); %pericentre radius of Earth
nEa  = sqrt(mu_Sun/kEa1(1)^3);  %medium angular velocity of the Earth
r_SUN=norm(rEa1); %[km] distance of Earth wrt Sun in 2/12/1995 
rSOI1=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth

% Position of L1 at departure date:
f = @(b) 1+(mu_Earth/mu_Sun)/((1-b)*b^2)-1/(1-b)^3 ;
b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)

%vEa1=vEa1*0.9;
rLag1=rEa1*(1-b);%[km] L1 position wrt to Sun on 2/12/1995
vLag1=vEa1*(1-b);%[km/s] L1 velocity wrt to Sun on 2/12/1995
distance_ratio1=norm(b*rEa1)/rSOI1; %ratio between distance of Earth from L1 and Earth SOI radius

[time, r, v] = extract_horizons('horizons_results.txt');

date_starting=jd2date(time(1));
date_ending=jd2date(time(end));


%FROM TIMELINE, THE MANEUVERS WERE MADE IN:

date1=[1995, 12, 3, 12, 0, 0];   %first orbit correction maneuvre
date2=[1995, 12, 22, 0, 0, 0];   %second midcourese correction maneuvre
date3=[1995, 12, 24, 8, 0, 0];   %nominal transfer orbit mode #
date4=[1996, 3, 26, 0, 0, 0];    %Halo Orbit Insertion maneuver
date5=[1996, 3, 30, 8, 0, 0];    %Halo Orbit

%FROM https://soho.nascom.nasa.gov/soc/soho_events/SOHO-Spacecraft-Events.pdf THE MANEUVERS WERE MADE IN:

DATE1=[1995, 12, 3, 0, 0, 0]; 
DATE2=[1996, 1, 5, 5, 40, 0]; % MCC2 X-1 Burn
DATE3=[1996, 2, 14, 17, 0, 0]; % Halo Orbit Insertion Manoeuvre (HOI)
DATE4=[1996, 3, 20, 22, 45, 0]; % Halo Orbit Insertion Manoeuvre (HOI) - Trim

%FROM THEORY:
DA1=[1995, 12, 2, 8, 8, 1]; %DEPARTURE DATE
TI1=date2jd(DA1); %DEPARTURE DATE IN JULIAN DAYS

%DeltaV to escape Earth from Cape Canaveral
ve0=sqrt(2*mu_Earth/r_Earth);%[km/s] overall deltaV necessary to escape Earth Gravity
v0=cos(lCCAS)*r_Earth*w_E;%[km/s] Cape Canaveral tangential velocity due to Earth spin rate
dv0=ve0-v0;%[km/s] deltaV expense to escape Earth Gravity
theta_SOI = acos((r_Earth/rSOI1)-1);   %HYPOTHESIS OF PARABOLIC DEPARTURE FROM EARTH
ToF_SOI_parabolic=ToF(Inf, 1, pi/2, theta_SOI, mu_Earth, r_Earth);
ToF_SOI_parabolic_days=ToF_SOI_parabolic/(24*3600); %TIME IN DAYS TO REACH THE SOI OF THE EARTH RADIUS IN DISTANCE

TI2=TI1+ToF_SOI_parabolic_days;
DA2=jd2date(TI2);

%set Three Body Problem
figure(1)
SAE=[-1.5*1e5*sin(23.45*pi/180), 0, 1.5*1e5*cos(23.45*pi/180)]'; %Earth spin axis
az_sae=2*pi*28.666667/365.25;
az_rr0=deg2rad(-150.93);                     %-158    %-150.93
R_sae=[cos(az_sae), sin(az_sae), 0;
      -sin(az_sae), cos(az_sae), 0;
       0,           0,           1];
SAE=R_sae*SAE;
plot3([0, SAE(1)], [0, SAE(2)], [0, SAE(3)], 'b', DisplayName='Earth spin axis');
hold on
rr0=r_Earth*[cos(az_rr0), sin(az_rr0), 0];
rr0(3)=-(SAE(1)*rr0(1)+SAE(2)*rr0(2))/SAE(3);
rr0=r_Earth*rr0/norm(rr0);
rr0=cos(lCCAS)*rr0+sin(lCCAS)*SAE'*r_Earth/norm(SAE);
rr0=r_Earth*rr0/norm(rr0);
vvt=cross(SAE, rr0)/norm(cross(SAE, rr0)); %tangential velocity direction from earth rotation
ar_vv0=deg2rad(13.02); %angle of rotation on the radial axis from the tangential direction                                              %8.8  %13.02
vv0_versor=vvt*cos(ar_vv0)+sin(ar_vv0)*cross(rr0, vvt)/norm(cross(rr0, vvt)); %initial velocity direction projected on the earth surface
vv0_plot=1.5*1e5*vv0_versor;
plot3([0, vv0_plot(1)], [0, vv0_plot(2)], [0, vv0_plot(3)], 'k', DisplayName='Initial Velocity direction');
hold on
vv0=ve0*vv0_versor; %actual initial velocity
vCCAS=v0*vvt; %cape canaveral velocity
t1=14; %ToF_SOI_parabolic_days;
dt1=50;
% rr0=-r_Earth*[0, -1, 0];
% alpha=321*pi/180;
% vv0=[0, cos(alpha), sin(alpha)];

[RRT, VVT, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on

%TRY LISSAJOUS ORBIT:
rrL1(1)=b*norm(rEa1)-1*Ax;  %1 %0.98
rrL1(2)=0;
rrL1(3)=-1*Az;  %1 %0.98
dt2=250;
T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=Av_distance/T_SOHO_s;
vvL1=[0, 1, 0]*(Ay/sqrt(Ax^2+Az^2))*1.12559043945*Av_velocity;   %0.98175  %0.9745
k=0.1; %[#] number of periods

[RRL1, VVL1, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rrL1, vvL1, rEa1, vEa1, k*T_SOHO, dt2);

co=3713;
coo=co+270;
date_analyzed=jd2date(time(coo));
%[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, r(:, coo)', v(:, coo)', rEa1, vEa1, k*T_SOHO, dt);

% find distance between the two trajectories (transfer traj. and lissajous traj.)
[m, ~]=size(RRT);
[n, ~]=size(RRL1);

p=25;
q=125;

rrdiff=zeros(floor(m/p), floor(n/q));
for i=1:floor(m/p)
    for j=1:floor(n/q)
        rrdiff(i, j)=norm(RRT(i*p, :)-RRL1(j*q, :));
    end
end
min(min(rrdiff))
[I, J]=find(rrdiff==min(min(rrdiff)));
I=I*p;
J=J*q;

RRT_INJ=RRT(I, :);
RRL1_INJ=RRL1(J, :);
VVT_INJ=VVT(I, :);
VVL1_INJ=VVL1(J, :);

dv1=norm(vv0-vCCAS);
tof1=t1*I/m;

% SOHO SIMULATION:
close all;

[RRT, VVT, rr_INJ, vv_INJ1] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, tof1, dt1);
hold on

TI_INJ=TI1+tof1;
DA_INJ=jd2date(TI_INJ);
%coordinates of Earth at injection day:
[kEa2,   ~] = uplanet(date2mjd2000(DA_INJ), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995

k2=2.5;
vv_INJ2=0.96171267*vvL1;     %0.9721081104
[RR_INJ, VV_INJ, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr_INJ, vv_INJ2, rEa2, vEa2, k2*T_SOHO, dt2);


dvv2=vv_INJ2-vv_INJ1;
dv2=norm(dvv2);

dv=dv1+dv2
%%
clc;
[m, ~]=size(RRT);
h=zeros(3, m);
h_norm=zeros(1, m);

for i=1:m
    h(:, i)=cross(RRT(i, :), VVT(i, :));
    h_norm(i)=norm(h(:, i));
end
figure(3)
plot(0:dt1:tof1, h_norm);
%% SOHO MISSION SIMULATION WITH THREE BODY PROBLEM SOLVING


clc; clear;
close all;
addpath(genpath('Default_Functions\'));
addpath(genpath('timeConversion\'));

G=6.6743e-20; %[km^3/(kg*s)]

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

% Ax = 200000;
% Ay = 650000;
% Az = 200000;


% Coordinates of Earth at departure of the mission:
[kEa1,   ~] = uplanet(date2mjd2000([1995, 12, 2, 8, 8, 1]), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa1, vEa1] = kep2car(kEa1(1), kEa1(2), kEa1(3), kEa1(4), kEa1(5), kEa1(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
rpEa = kEa1(1)*(1-kEa1(2)); %pericentre radius of Earth
nEa  = sqrt(mu_Sun/kEa1(1)^3);  %medium angular velocity of the Earth
r_SUN=norm(rEa1); %[km] distance of Earth wrt Sun in 2/12/1995 
rSOI1=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth

 

% Position of L1 at departure date:
f = @(b) 1+(mu_Earth/mu_Sun)/((1-b)*b^2)-1/(1-b)^3 ;
b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)

%vEa1=vEa1*0.9;
rLag1=rEa1*(1-b);%[km] L1 position wrt to Sun on 2/12/1995
vLag1=vEa1*(1-b);%[km/s] L1 velocity wrt to Sun on 2/12/1995
distance_ratio1=norm(b*rEa1)/rSOI1; %ratio between distance of Earth from L1 and Earth SOI radius

[time, r, v] = extract_horizons('horizons_results.txt');

date_starting=jd2date(time(1));
date_ending=jd2date(time(end));


%FROM TIMELINE, THE MANEUVERS WERE MADE IN:

date1=[1995, 12, 3, 12, 0, 0];   %first orbit correction maneuvre
date2=[1995, 12, 22, 0, 0, 0];   %second midcourese correction maneuvre
date3=[1995, 12, 24, 8, 0, 0];   %nominal transfer orbit mode #
date4=[1996, 3, 26, 0, 0, 0];    %Halo Orbit Insertion maneuver
date5=[1996, 3, 30, 8, 0, 0];    %Halo Orbit

%FROM https://soho.nascom.nasa.gov/soc/soho_events/SOHO-Spacecraft-Events.pdf THE MANEUVERS WERE MADE IN:

DATE1=[1995, 12, 3, 0, 0, 0]; 
DATE2=[1996, 1, 5, 5, 40, 0]; % MCC2 X-1 Burn
DATE3=[1996, 2, 14, 17, 0, 0]; % Halo Orbit Insertion Manoeuvre (HOI)
DATE4=[1996, 3, 20, 22, 45, 0]; % Halo Orbit Insertion Manoeuvre (HOI) - Trim

%FROM THEORY:
DA1=[1995, 12, 2, 8, 8, 1]; %DEPARTURE DATE
TI1=date2jd(DA1); %DEPARTURE DATE IN JULIAN DAYS

%DeltaV to escape Earth from the parking orbit
r0=r_Earth+180; %[km] radius of the initial parking orbit
ve0=sqrt(2*mu_Earth/r0);%[km/s] overall velocity necessary to escape Earth Gravity
v0=sqrt(mu_Earth/r0);%[km/s] velocity in the parking orbit

%set Three Body Problem
figure(1)
eps=deg2rad(23.45); %[rad] inclination of the Earth spin axis wrt ecliptic
i=deg2rad(28.8); %[rad] inclination of the initial parking orbit
az_sae=deg2rad(-28.666667*360/365.25);    %[rad] rotation of the earth spin axis wrt the z axis of the inertial frame
az_rr0=deg2rad(-159);                     %-158    %-150.93
ar_vv0=deg2rad(0); %angle of rotation OF THE VELOCITY on the radial axis from the tangential direction                                              %8.8  %13.02
R_sae=[cos(az_sae), sin(az_sae), 0;
      -sin(az_sae), cos(az_sae), 0;
       0,           0,           1];
SAE=[-1.5*1e5*sin(eps), 0, 1.5*1e5*cos(eps)]'; %Earth spin axis
SAE=R_sae*SAE;
plot3([0, SAE(1)], [0, SAE(2)], [0, SAE(3)], 'b', DisplayName='Earth spin axis');
hold on
GAM=1.5*1e5*R_sae*[0; 1; 0];
plot3([0, GAM(1)], [0, GAM(2)], [0, GAM(3)], 'r', DisplayName='Gamma line');
hold on
H0=[-1.5*1e5*sin(eps-i), 0, 1.5*1e5*cos(eps-i)]'; %Earth spin axis
H0=R_sae*H0;
plot3([0, H0(1)], [0, H0(2)], [0, H0(3)], 'm', DisplayName='Initial momentum');
hold on
rr0=r0*[cos(az_rr0), sin(az_rr0), 0];  
rr0(3)=-(H0(1)*rr0(1)+H0(2)*rr0(2))/SAE(3);
rr0=r0*rr0/norm(rr0);                  %initial position on the parking orbit
vvt=cross(H0, rr0)/norm(cross(H0, rr0)); %tangential velocity direction from INITIAL PARKING ORBIT
vv0_versor=vvt*cos(ar_vv0)+sin(ar_vv0)*cross(rr0, vvt)/norm(cross(rr0, vvt)); %initial velocity direction FOR ESCAPE
vv0_plot=1.5*1e5*vv0_versor; 
plot3([0, vv0_plot(1)], [0, vv0_plot(2)], [0, vv0_plot(3)], 'k', DisplayName='Initial Velocity direction');
hold on
vv0=ve0*vv0_versor; %actual initial velocity FOR ESCAPE
vPO=v0*vvt; %PARKING ORBIT VELOCITY
t1=14; %ToF_SOI_parabolic_days;
dt1=250;


[RRT, VVT, rr_b, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on
close all;
[m, ~]=size(rr_b);
I=find(abs(rr_b(:, 3))>Az & rr_b(:, 1)>11*1e5);

J=find(abs(rr_b(:, 2))==min(abs(rr_b(I(1):I(end), 2))));


TOF=t1*J/m;

[RRT, VVT, rr_b, VVB] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, TOF, dt1);

% Position of Earth at arrival date:
TIH=TI1+TOF;
DAH=jd2date(TIH);
%coordinates of Earth at injection day:
[kEa2,   ~] = uplanet(date2mjd2000(DAH), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=Av_distance/T_SOHO_s;
vvH=[0, 1, 0]*(Ay/sqrt(Ax^2+Az^2))*0.94539615*Av_velocity;
k=2; %[#] number of periods
dt2=750;

[RRH, VVH, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr_b(end, :), vvH, rEa2, vEa2, k*T_SOHO, dt2);
%%
%TRY LISSAJOUS ORBIT:
rrL1(1)=b*norm(rEa1)-1*Ax;  %1 %0.98
rrL1(2)=0;
rrL1(3)=-1*Az;  %1 %0.98
dt2=750;
T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=Av_distance/T_SOHO_s;
vvL1=[0, 1, 0]*(Ay/sqrt(Ax^2+Az^2))*1.12559043945*Av_velocity;   %0.98175  %0.9745
k=1; %[#] number of periods

[RRL1, VVL1, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rrL1, vvL1, rEa1, vEa1, k*T_SOHO, dt2);

co=3713;
coo=co+270;
date_analyzed=jd2date(time(coo));
%[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, r(:, coo)', v(:, coo)', rEa1, vEa1, k*T_SOHO, dt);

% find distance between the two trajectories (transfer traj. and lissajous traj.)
[m, ~]=size(RRT);
[n, ~]=size(RRL1);

p=5;
q=25;

rrdiff=zeros(floor(m/p), floor(n/q));
for i=1:floor(m/p)
    for j=1:floor(n/q)
        rrdiff(i, j)=norm(RRT(i*p, :)-RRL1(j*q, :));
    end
end
min(min(rrdiff))
[I, J]=find(rrdiff==min(min(rrdiff)));
I=I*p;
J=J*q;

RRT_INJ=RRT(I, :);
RRL1_INJ=RRL1(J, :);
VVT_INJ=VVT(I, :);
VVL1_INJ=VVL1(J, :);

dv1=norm(vv0-vPO);
tof1=t1*I/m;

%% SOHO SIMULATION:
close all;

[RRT, VVT, rr_INJ, vv_INJ1] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, tof1, dt1);
hold on

TI_INJ=TI1+tof1;
DA_INJ=jd2date(TI_INJ);
%coordinates of Earth at injection day:
[kEa2,   ~] = uplanet(date2mjd2000(DA_INJ), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995

k2=3;
vv_INJ2=0.9721081104*vvL1;     %0.9721081104
[RR_INJ, VV_INJ, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr_INJ, vv_INJ2, rEa2, vEa2, k2*T_SOHO, dt2);


dvv2=vv_INJ2-vv_INJ1;
dv2=norm(dvv2);

dv=dv1+dv2


%% SOHO MISSION SIMULATION WITH THREE BODY PROBLEM SOLVING


clc; clear;
close all;
addpath(genpath('Default_Functions\'));
addpath(genpath('timeConversion\'));

G=6.6743e-20; %[km^3/(kg*s)]

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

% Ax = 200000;
% Ay = 650000;
% Az = 200000;


% Coordinates of Earth at departure of the mission:
[kEa1,   ~] = uplanet(date2mjd2000([1995, 12, 2, 8, 8, 1]), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa1, vEa1] = kep2car(kEa1(1), kEa1(2), kEa1(3), kEa1(4), kEa1(5), kEa1(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
rpEa = kEa1(1)*(1-kEa1(2)); %pericentre radius of Earth
nEa  = sqrt(mu_Sun/kEa1(1)^3);  %medium angular velocity of the Earth
r_SUN=norm(rEa1); %[km] distance of Earth wrt Sun in 2/12/1995 
rSOI1=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth

 

% Position of L1 at departure date:
f = @(b) 1+(mu_Earth/mu_Sun)/((1-b)*b^2)-1/(1-b)^3 ;
b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)

%vEa1=vEa1*0.9;
rLag1=rEa1*(1-b);%[km] L1 position wrt to Sun on 2/12/1995
vLag1=vEa1*(1-b);%[km/s] L1 velocity wrt to Sun on 2/12/1995
distance_ratio1=norm(b*rEa1)/rSOI1; %ratio between distance of Earth from L1 and Earth SOI radius

[time, r, v] = extract_horizons('horizons_results.txt');

date_starting=jd2date(time(1));
date_ending=jd2date(time(end));


%FROM TIMELINE, THE MANEUVERS WERE MADE IN:

date1=[1995, 12, 3, 12, 0, 0];   %first orbit correction maneuvre
date2=[1995, 12, 22, 0, 0, 0];   %second midcourese correction maneuvre
date3=[1995, 12, 24, 8, 0, 0];   %nominal transfer orbit mode #
date4=[1996, 3, 26, 0, 0, 0];    %Halo Orbit Insertion maneuver
date5=[1996, 3, 30, 8, 0, 0];    %Halo Orbit

%FROM https://soho.nascom.nasa.gov/soc/soho_events/SOHO-Spacecraft-Events.pdf THE MANEUVERS WERE MADE IN:

DATE1=[1995, 12, 3, 0, 0, 0]; 
DATE2=[1996, 1, 5, 5, 40, 0]; % MCC2 X-1 Burn
DATE3=[1996, 2, 14, 17, 0, 0]; % Halo Orbit Insertion Manoeuvre (HOI)
DATE4=[1996, 3, 20, 22, 45, 0]; % Halo Orbit Insertion Manoeuvre (HOI) - Trim

%FROM THEORY:
DA1=[1995, 12, 2, 8, 8, 1]; %DEPARTURE DATE
TI1=date2jd(DA1); %DEPARTURE DATE IN JULIAN DAYS

%DeltaV to escape Earth from the parking orbit
r0=r_Earth+180; %[km] radius of the initial parking orbit
ve0=sqrt(2*mu_Earth/r0);%[km/s] overall velocity necessary to escape Earth Gravity
v0=sqrt(mu_Earth/r0);%[km/s] velocity in the parking orbit

%set Three Body Problem
figure(1)
eps=deg2rad(23.45); %[rad] inclination of the Earth spin axis wrt ecliptic
i=deg2rad(28.8); %[rad] inclination of the initial parking orbit
az_sae=deg2rad(-28.666667*360/365.25);    %[rad] rotation of the earth spin axis wrt the z axis of the inertial frame
az_rr0=deg2rad(-159);                     %-158    %-150.93
ar_vv0=deg2rad(0); %angle of rotation OF THE VELOCITY on the radial axis from the tangential direction                                              %8.8  %13.02
R_sae=[cos(az_sae), sin(az_sae), 0;
      -sin(az_sae), cos(az_sae), 0;
       0,           0,           1];
SAE=[-1.5*1e5*sin(eps), 0, 1.5*1e5*cos(eps)]'; %Earth spin axis
SAE=R_sae*SAE;
plot3([0, SAE(1)], [0, SAE(2)], [0, SAE(3)], 'b', DisplayName='Earth spin axis', HandleVisibility='off');
hold on
GAM=1.5*1e5*R_sae*[0; 1; 0];
plot3([0, GAM(1)], [0, GAM(2)], [0, GAM(3)], 'r', DisplayName='Gamma line', HandleVisibility='off');
hold on
H0=[-1.5*1e5*sin(eps-i), 0, 1.5*1e5*cos(eps-i)]'; %Earth spin axis
H0=R_sae*H0;
plot3([0, H0(1)], [0, H0(2)], [0, H0(3)], 'm', DisplayName='Initial momentum', HandleVisibility='off');
hold on
rr0=r0*[cos(az_rr0), sin(az_rr0), 0];  
rr0(3)=-(H0(1)*rr0(1)+H0(2)*rr0(2))/SAE(3);
rr0=r0*rr0/norm(rr0);                  %initial position on the parking orbit
vvt=cross(H0, rr0)/norm(cross(H0, rr0)); %tangential velocity direction from INITIAL PARKING ORBIT
vv0_versor=vvt*cos(ar_vv0)+sin(ar_vv0)*cross(rr0, vvt)/norm(cross(rr0, vvt)); %initial velocity direction FOR ESCAPE
vv0_plot=1.5*1e5*vv0_versor; 
plot3([0, vv0_plot(1)], [0, vv0_plot(2)], [0, vv0_plot(3)], 'k', DisplayName='Initial Velocity direction', HandleVisibility='off');
hold on
vv0=ve0*vv0_versor; %actual initial velocity FOR ESCAPE
vPO=v0*vvt; %PARKING ORBIT VELOCITY
t1=14; %ToF_SOI_parabolic_days;
dt1=250;


[RRT, VVT, rr_b, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on
close all;
[m, ~]=size(rr_b);
I=find(abs(rr_b(:, 3))>Az & rr_b(:, 1)>11*1e5);

J=find(abs(rr_b(:, 2))==min(abs(rr_b(I(1):I(end), 2))));


TOF=t1*J/m;

[RRT, VVT, rr_b, VVB] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, TOF, dt1);

% Position of Earth at arrival date:
TIH=TI1+TOF;
DAH=jd2date(TIH);
%coordinates of Earth at injection day:
[kEa2,   ~] = uplanet(date2mjd2000(DAH), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=Av_distance/T_SOHO_s;
vvH=[0, 1, 0]*(Ay/sqrt(Ax^2+Az^2))*0.94539615*Av_velocity;
k=2; %[#] number of periods
dt2=750;

[RRH, VVH, ~, ~] = TBP2_2(mu_Sun/G, mu_Earth/G, rr_b(end, :), vvH, rEa2, vEa2, k*T_SOHO, dt2);

dvTTI=ve0-v0; %dv for insertion in quasi-parabolic transfer orbit
dvHOI=norm(vvH-VVB); %dv for insertion in halo lissajous orbit