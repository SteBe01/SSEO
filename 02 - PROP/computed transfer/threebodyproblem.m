clc;
clear;
d=384000000;

M=5.972e+18;
m=7.348e+22;
t=2e9;
dt=300;
x=6.25;
y=0.5*x;
rr0=[-1e7; 0;    0];
vv0=[   0; x;   x];

TerraLuna(M, M, d, rr0, vv0, t, dt, 6378000)

%%
clc;
clear;
d=384000000;

M=5.972e+18;
m=7.348e+22;
t=2e9;
dt=300;
x=6.3;
y=0.5*x;
rr0=[-1e7; 0;    0];
vv0=[   0; x;   x];

TerraLuna(M, M, d, rr0, vv0, t, dt, 6378000)

%%
clc;
clear;
d=384000000;

M=5.972e+18;
m=7.348e+22;
t=2e9;
dt=300;
x=6.1;
y=0.5*x;
rr0=[-1e7; 0;    0];
vv0=[   0; x;   x];

TerraLuna(M, M, d, rr0, vv0, t, dt, 6378000)

%%
clc;
clear;
d=384000000;

M=5.972e+18;
m=7.348e+22;
t=1e8;
dt=1000;
x=0.01;
y=5*x;
rr0=[ d/2-1; 0;    0];
vv0=[  x; -y;   y];

TerraLuna(M, M, d, rr0, vv0, t, dt, 6378000)

%%
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e4;
dt=1;

rr0=[2310687.3; -8071383.9; -2384194.2];
vv0=[   6124.0;     1924.0;    -1176.0];

TerraLuna(M, m, d, rr0, vv0, t, dt)

%%
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1.4e6;
dt=1;

rr0=[2310687.3; -8071383.9; -2384194.2];
vv0=1.47*[   6124.0;     1924.0;    -1176.0];

TreCorpi(M, m, d, rr0, vv0, t, dt, 'L')

%%
clc;
clear;
d=38440000;

M=5.972e+24;
m=7.348e+22;
t=5.1e5;
dt=1;

rr0=[2310687.3; -8071383.9; -2384194.2];
vv0=1.06*[   6124.0;     1924.0;    -1176.0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000)

%% orbita stabile
clc;
clear;
d=38440000;

M=5.972e+24;
m=7.348e+22;
t=(7/25)*1e5;
dt=0.1;

rr0=1e7*[-1; 0; 0];
vv0=1e4*[0; 1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000)

%% quasi orbita a 8
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e6;
dt=1;

rr0=[d/2; 0; 0];
vv0=0.225e4*[1; 1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% fionda gravitazionale
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=(5/25)*1e6;
dt=1;

rr0=1e7*[-1; 0; 0];
vv0=0.9e4*[0; 1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000)

%% 
clc;
clear;
d=38440000;

M=5.972e+24;
m=7.348e+22;
t=1e5;
dt=1;


th=50;

th1=th*pi/180;

v=2;

rr0=[d/2; 0; 0];
vv0=v*[cos(th1); sin(th1); 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000)

%% 
clc;
clear;
d=3844000000;

M=5.972e+24;
m=7.348e+22;
t=1e5;
dt=0.1;

th=90;

th1=th*pi/180;

v=-6e3;

rr0=1e7*[-1; 0; 0];
vv0=v*[cos(th1); sin(th1); 0];

TreCorpi(M, m, d, rr0, vv0, t, dt)

%% orbita finale IAMS

clc;
clear;

a2=13270000;
e2=0.3469;
i2=1.4580;
OM2=2.5710;
om2=1.8860;
theta2=1.9710;
[rr2, vv2]=parorb2rv_p(a2, e2, i2, OM2, om2, theta2);

d=3844000000;

M=5.972e+24;
m=7.348e+22;
t=1e5;
dt=0.1;
TreCorpi(M, m, d, rr2, vv2, t, dt)

%% quasi orbita a 8
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e6;
dt=1;

rr0=[d/2; 0; 0];
vv0=0.2e4*[1; 1.2; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% 
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e6;
dt=1;

rr0=[d/2; 0; 0];
vv0=0.1e3*[1; 3; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% 
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e5;
dt=0.1;
th=30;

th1=th*pi/180;

v=4e3;

rr0=[d/2; 0; 0];
vv0=v*[cos(th1); sin(th1); 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% quasi traiettoria a 8
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=4e5;
dt=0.1;
th=24;

th1=th*pi/180;

v=4e3;

rr0=[d/2; 0; 0];
vv0=v*[cos(th1); sin(th1); 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% 
clc;
clear;
d=384400000;

M=5.972e+24;
m=7.348e+22;
t=1e6;
dt=0.1;
th=80;

th1=th*pi/180;

v=5e2;

rr0=[d/2; -d/20; 0];
vv0=v*[cos(th1); sin(th1); 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% orbita curiosa 
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=1e6;
dt=0.1;
x=d/2;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr;
vv0=v*[0; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% orbita curiosa 
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=2e6;
dt=1;
x=d/2;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr+100
vv0=v*[0; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% orbita curiosa 
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=5e5;
dt=0.1;
x=d/4;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr*2;
vv0=v*[0; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% fionda gravitazionale 
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=5e5;
dt=1;
x=d/4;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr*2.2;
vv0=v*[0.1; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% orbita curiosa
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=4e6;
dt=1;
x=d*0.53;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr;
vv0=v*[0; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% 
clc;
clear;
d=384400000;
G=6.6743e-11;
M=5.972e+24;

t=1e5;
dt=0.1;
x=d*0.826;
rr0=[-x; 0; 0];

om=sqrt((G*2*M)/(d^3));
vr=om*(x+d/2);
vo=sqrt(G*M/x);
v=vr;
vv0=v*[0; -1; 0];

TreCorpi(M, M, d, rr0, vv0, t, dt, 6378000);

%% 20/03/2024 TBP fionda gravitazionale

clc; clear;
close all;
addpath(genpath('Functions'));

G=6.6743e-20; %[km^3/(kg*s)]

mu_Sun = astroConstants(4); %[km^3/s^2] Sun gravity constant
mu_Earth = astroConstants(13)*100000; %[km^3/s^2] Earth gravity constant
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

t=200;
dt=500;

% rr0=rEa1*0.05;
% rr0(3)=1000000;
rr0=-rEa1+norm(rEa1)*[0, 0, 1]';
vv0=vEa1*1e-10;

[RRT, VVT] = TBP(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1*1e0, vEa1*1e0, t, dt);


%% 22/03/2024 TBP - L1 point verification - instability for excess of velocity

clc; clear;
close all;
addpath(genpath('Functions'));

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

%set Three Body Problem
t=800;
dt=500;

rr0=-rEa1+rLag1;
vv0=-vEa1+vLag1+0.0001*vEa1;

[RRT, VVT] = TBP(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t, dt);

%% 22/03/2024 TBP - L1 point verification - unstable

clc; clear;
close all;
addpath(genpath('Functions'));

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

%set Three Body Problem
t=900;
dt=500;

rr0=-rEa1+rLag1+0.000001*rEa1;
vv0=-vEa1+vLag1-0.000001*vEa1;

[RRT, VVT] = TBP(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t, dt);

%% 22/03/2024 TBP - L1 point verification

clc; clear;
close all;
addpath(genpath('Functions'));

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

%set Three Body Problem
t=900;
dt=500;

rr0=-rEa1+rLag1-0.000001*rEa1;
vv0=-vEa1+vLag1+0.00001*vEa1;

[RRT, VVT] = TBP(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t, dt);
%% 22/03/2024 TBP - L2 point verification

clc; clear;
close all;
addpath(genpath('Functions'));

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
rLag2=rEa1*(1+b);%[km] L2 position wrt to Sun on 2/12/1995
vLag2=vEa1*(1+b);%[km/s] L2 velocity wrt to Sun on 2/12/1995
distance_ratio1=norm(b*rEa1)/rSOI1; %ratio between distance of Earth from L1 and Earth SOI radius

%set Three Body Problem
t=400;
dt=500;

rr0=-rEa1+rLag2;
vv0=-vEa1+vLag2+0.0001*vEa1;

[RRT, VVT] = TBP(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t, dt);


%% 24/03/2024 TBP - SOHO mission path simulation - halo orbit stable for 3 periods

clc; clear;
close all;
addpath(genpath('Functions'));

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
t1=18; %ToF_SOI_parabolic_days;
dt1=250;

rr0=-r_Earth*[1, 0, 0];
alpha=321*pi/180;
vv0=-ve0*[0, cos(alpha), sin(alpha)];

[RRT, VVT] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on



%TRY LISSAJOUS ORBIT:
rrL1(1)=b*norm(rEa1)-1*Ax;  %1 %0.98
rrL1(2)=0;
rrL1(3)=-1*Az;  %1 %0.98
dt2=500;
T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=Av_distance/T_SOHO_s;
vvL1=[0, 1, 0]*(Ay/sqrt(Ax^2+Az^2))*1.12559043945*Av_velocity;   %0.98175  %0.9745
k=3; %[#] number of periods

[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, rrL1, vvL1, rEa1, vEa1, k*T_SOHO, dt2);

co=3713;
coo=co+270;
date_analyzed=jd2date(time(coo));
%[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, r(:, coo)', v(:, coo)', rEa1, vEa1, k*T_SOHO, dt);


%% find distance between the two trajectories (transfer traj. and lissajous traj.)
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

%% 24/03/2024 TBP - try to find a stable orbit around L1

clc; clear;
close all;
addpath(genpath('Functions'));

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
t1=18; %ToF_SOI_parabolic_days;
dt=250;

rr0=-r_Earth*[1, 0, 0];
alpha=141*pi/180;
vv0=-ve0*[0, cos(alpha), sin(alpha)];

%[RRT, VVT] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt);
%hold on





T_SOHO=178; %[days] period of SOHO S/C around L1
T_SOHO_s=T_SOHO*3600*24; %[s] period of SOHO S/C around L1


co=3713;
coo=co+270;
date_analyzed=jd2date(time(coo));
%[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, r(:, coo)', v(:, coo)', rEa1, vEa1, k*T_SOHO, dt);


%Set problem constants:
omL1=2*pi/(T_SOHO_s);        %angular velocity in lissajous orbit around L1 
omEa=norm(vEa1)/norm(rEa1);  %angular velocity of the system on 2/12/1995
ua=norm(rEa1);
f1 = @(x1, x2) mu_Earth/((x1^2+x2^2)^1.5) + mu_Sun/((x1^2+(1-x2)^2)^1.5) - omL1^2 * ua^3;
f2 = @(x1, x2) (omEa^2)*(1-x2)*ua^2 + x2*mu_Earth/(ua *(x1^2+x2^2)^1.5) - (1-x2)*mu_Sun/(ua *(x1^2+(1-x2)^2)^1.5);

vec=[b-0.0045:1e-6:b+0.00001];
k=1;
for x2=b %0.009:1e-6:0.01
    f1_ = @(x) f1(x, x2);
    f2_ = @(x) f2(x, x2);

    x1_field=[0, 0.5];
    X1_1(k)=fzero(f1_, x1_field);
    X1_2(k)=fzero(f2_, 0.01);
    X1_D(k)=abs(X1_2(k)-X1_1(k));
    k=k+1;
end
figure(3)
plot(vec, X1_D);

% 398600/(x^2+y^2)^1,5+132712440018/(x^2+(1-y)^2)^1,5- 147493744^3*(4,08550854867587*10^-7)^2
% y*398600/(x^2+y^2)^1,5+(1-y)*132712440018/(x^2+(1-y)^2)^1,5- 147493744^3*(2,04799272949247*10^-7)^2


%TRY LISSAJOUS ORBIT:
rrL1(1)=b*ua;
rrL1(2)=0;
rrL1(3)=-X1_1*ua;
Av_distance=pi*(Ay+sqrt(Ax^2+Az^2)); %[km] average distance crossed by SOHO S/C in his orbit every period around L1
Av_velocity=X1_1*ua*omL1;
vvL1=[0, 1, 0]*Av_velocity;   
k=2; %[#] number of periods

[RRL1, VVL1] = TBP2(mu_Sun/G, mu_Earth/G, rrL1, vvL1, rEa1, vEa1, k*T_SOHO, dt);


%% 23/03/2024 TBP - SOHO mission path simulation

clc; clear;
close all;
addpath(genpath('Functions'));
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
ar=deg2rad(13.02); %angle of rotation on the radial axis from the tangential direction                                              %8.8  %13.02
vv0_versor=vvt*cos(ar)+sin(ar)*cross(rr0, vvt)/norm(cross(rr0, vvt)); %initial velocity direction projected on the earth surface
vv0_plot=1.5*1e5*vv0_versor;
plot3([0, vv0_plot(1)], [0, vv0_plot(2)], [0, vv0_plot(3)], 'k', DisplayName='Initial Velocity direction');
hold on
vv0=ve0*vv0_versor; %actual initial velocity
vCCAS=v0*vvt; %cape canaveral velocity
t1=14; %ToF_SOI_parabolic_days;
dt1=250;
% rr0=-r_Earth*[0, -1, 0];
% alpha=321*pi/180;
% vv0=[0, cos(alpha), sin(alpha)];

[RRT, VVT, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on

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

dv1=norm(vv0-vCCAS);
% dvv2=VVL1_INJ-VVT_INJ;
% dv2=norm(VVL1_INJ-VVT_INJ);


tof1=t1*I/m

% dv=dv1+dv2

%% SOHO SIMULATION
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



%% 26/03/2024 TBP - Simple Earth to Earth-Sun L1 transfer

clc; clear;
close all;
addpath(genpath('Functions'));

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

TOF=74; %[days]

%FROM THEORY:
DA1=[1995, 12, 2, 8, 8, 1]; %DEPARTURE DATE
TI1=date2jd(DA1); %DEPARTURE DATE IN JULIAN DAYS
TI2=TI1+TOF;
DA2=jd2date(TI2);

% Coordinates of Earth at departure of the mission:
[kEa1,   ~] = uplanet(date2mjd2000(DA1), 3) ;  %get the keplearian coordinates of Earth on Launch date and hour (2/12/1995)
[rEa1, vEa1] = kep2car(kEa1(1), kEa1(2), kEa1(3), kEa1(4), kEa1(5), kEa1(6), mu_Sun); %cartesian coordinates of the Earth at 2/12/1995
[kEa2,   ~] = uplanet(date2mjd2000(DA2), 3) ;  %get the keplearian coordinates of Earth on Arrival date and hour
[rEa2, vEa2] = kep2car(kEa2(1), kEa2(2), kEa2(3), kEa2(4), kEa2(5), kEa2(6), mu_Sun); %cartesian coordinates of the Earth at Arrival
rpEa = kEa1(1)*(1-kEa1(2)); %pericentre radius of Earth
nEa  = sqrt(mu_Sun/kEa2(1)^3);  %medium angular velocity of the Earth
r_SUN=norm(rEa2); %[km] distance of Earth wrt Sun at arrival 
rSOI2=r_SUN*(mu_Earth/mu_Sun)^(2/5); %[km] radius of the Sphere Of Influence of the Earth

% Position of L1 at departure and arrival date:
f = @(b) 1+(mu_Earth/mu_Sun)/((1-b)*b^2)-1/(1-b)^3 ;
b=fzero(f, [1e-6, 1-1e-6]); %ratio between distance Earth-L1 and distance Earth-Sun (solution of Lagrange problem)
rLag1=rEa1*(1-b);%[km] L1 position wrt to Sun on 2/12/1995
vLag1=vEa1*(1-b);%[km/s] L1 velocity wrt to Sun on 2/12/1995
rLag2=rEa2*(1-b);%[km] L1 position wrt to Sun at arrival
vLag2=vEa2*(1-b);%[km/s] L1 velocity wrt to Sun at arrival
distance_ratio1=norm(b*rEa2)/rSOI2; %ratio between distance of Earth from L1 and Earth SOI radius

% [time, r, v] = extract_horizons('horizons_results.txt');
% 
% date_starting=jd2date(time(1));
% date_ending=jd2date(time(end));
% 
%
% %FROM TIMELINE, THE MANEUVERS WERE MADE IN:
% 
% date1=[1995, 12, 3, 12, 0, 0];   %first orbit correction maneuvre
% date2=[1995, 12, 22, 0, 0, 0];   %second midcourese correction maneuvre
% date3=[1995, 12, 24, 8, 0, 0];   %nominal transfer orbit mode #
% date4=[1996, 3, 26, 0, 0, 0];    %Halo Orbit Insertion maneuver
% date5=[1996, 3, 30, 8, 0, 0];    %Halo Orbit
% 
% %FROM https://soho.nascom.nasa.gov/soc/soho_events/SOHO-Spacecraft-Events.pdf THE MANEUVERS WERE MADE IN:
% 
% DATE1=[1995, 12, 3, 0, 0, 0]; 
% DATE2=[1996, 1, 5, 5, 40, 0]; % MCC2 X-1 Burn
% DATE3=[1996, 2, 14, 17, 0, 0]; % Halo Orbit Insertion Manoeuvre (HOI)
% DATE4=[1996, 3, 20, 22, 45, 0]; % Halo Orbit Insertion Manoeuvre (HOI) - Trim


%DeltaV to escape Earth from Cape Canaveral
ve0=sqrt(2*mu_Earth/r_Earth);%[km/s] overall deltaV necessary to escape Earth Gravity
v0=cos(lCCAS)*r_Earth*w_E;%[km/s] Cape Canaveral tangential velocity due to Earth spin rate
dv0=ve0-v0;%[km/s] deltaV expense to escape Earth Gravity
theta_SOI = acos((r_Earth/rSOI2)-1);   %HYPOTHESIS OF PARABOLIC DEPARTURE FROM EARTH
ToF_SOI_parabolic=ToF(Inf, 1, pi/2, theta_SOI, mu_Earth, r_Earth);
ToF_SOI_parabolic_days=ToF_SOI_parabolic/(24*3600); %TIME IN DAYS TO REACH THE SOI OF THE EARTH RADIUS IN DISTANCE

TI2=TI1+ToF_SOI_parabolic_days;
DA2=jd2date(TI2);

%set Three Body Problem
figure(1)
% SAE=[-1.5*1e5*sin(23.45*pi/180), 0, 1.5*1e5*cos(23.45*pi/180)]'; %Earth spin axis
% az_sae=2*pi*28.666667/365.25;
az_rr0=deg2rad(-157);                     %-158    %-150.93
% R_sae=[cos(az_sae), sin(az_sae), 0;
%       -sin(az_sae), cos(az_sae), 0;
%        0,           0,           1];
%SAE=R_sae*SAE;
% plot3([0, SAE(1)], [0, SAE(2)], [0, SAE(3)], 'b', DisplayName='Earth spin axis');
% hold on
rr0=r_Earth*[cos(az_rr0), sin(az_rr0), 0];
vvt=cross([0, 0, 1]', rr0)/norm(cross([0, 0, 1]', rr0)); %tangential velocity direction from earth rotation
% ar=deg2rad(13.02); %angle of rotation on the radial axis from the tangential direction                                              %8.8  %13.02
% vv0_versor=vvt*cos(ar)+sin(ar)*cross(rr0, vvt)/norm(cross(rr0, vvt)); %initial velocity direction projected on the earth surface
vv0_plot=1.5*1e5*vvt;
plot3([0, vv0_plot(1)], [0, vv0_plot(2)], [0, vv0_plot(3)], 'k', DisplayName='Initial Velocity direction');
hold on
vv0=ve0*vvt; %actual initial velocity
vCCAS=v0*vvt; %cape canaveral velocity
t1=1600; %ToF_SOI_parabolic_days;
dt1=250;
% rr0=-r_Earth*[0, -1, 0];
% alpha=321*pi/180;
% vv0=[0, cos(alpha), sin(alpha)];

[RRT, VVT, ~, ~] = TBP2(mu_Sun/G, mu_Earth/G, rr0, vv0, rEa1, vEa1, t1, dt1);
hold on