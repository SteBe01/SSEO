function [rr_t, vv_t] = TerraLuna(M, m, d, rr0, vv0, t, dt, Rl)
%TERRALUNA Summary of this function goes here
%   Detailed explanation goes here
G=6.6743e-11; %[m^3/(kg*s^2)]

if nargin<6

    dt=10;

end

X=[1; 0; 0];
Y=[0; 1; 0];
Z=[0; 0; 1];

DD=-X*(d*m/(M+m)); %posizione di M rispetto a O
dd=+X*(d*M/(M+m)); %posizione di m rispetto a O



%f=G*M/(d^2); %attrazione di m con M (accelerazione)
%om=sqrt(f/norm(dd)); 
om=sqrt((G*(M+m))/(d^3));%velocità angolare angolare del sistema omm=[0, 0, om];
omm=[0, 0, om];

tt=[0:dt:t];%vettore dei tempi
h=length(tt);

rr=zeros(3, h); %inizializzo rv
rr(:, 1)=rr0+DD; %posizione satellite rispetto a O al tempo 0

aa=zeros(3, h); %inizializzo accelerazioni satellite
vv=zeros(3, h); %inizializzo velocità (plurale) satellite

rrm=zeros(3, h);
rrM=zeros(3, h);
rt=zeros(3, h);
omR=zeros(1, h);
aom=zeros(3, h);
aam=zeros(3, h);
aaM=zeros(3, h);
dv=zeros(3, h);
%vr=zeros(3, h);

vvt=cross(omm, DD)'; %velocità di M rispetto a O

vv(:, 1)=vvt+vv0; %velocità del satellite rispetto a O

%vr(:, 1)=vv0; %velocità relativa del satellite rispetto a M

for k=1:h

    rrm(:, k)=rr(:, k)-dd; %posizione al passo k rispetto a m
    rrM(:, k)=rr(:, k)-DD; %posizione al passo k rispetto a M

    rt(:, k)=cross(Z, rr(:, k))/norm(rr(:, k)); %versore trasverso
    omR(k)=dot(rt(:, k), vv(:, k))/norm(rr(:, k)); %velocità angolare del satellite relativa rispetto alla terna fissa del sistema

    aom(:, k)=((om+omR(k))^2)*rr(:, k); %accelerazione centrifuga al passo k
    aam(:, k)=-G*m*rrm(:, k)/(norm(rrm(:, k))^3); %accelerazione grav per m al passo k
    aaM(:, k)=-G*M*rrM(:, k)/(norm(rrM(:, k))^3); %accelerazione grav per M al passo k

    aa(:, k)=aom(:, k)+aam(:, k)+aaM(:, k); %accelerazione complessiva al passo k

    dv(:, k)=aa(:, k)*dt; %delta di velocità tra passo k+1 e k

    vv(:, k+1)=vv(:, k)+dv(:, k); %velocità al passo k+1
    %vr(:, k+1)=vr(:, k)+dv(:, k); %velocità relativa al psso k+1

    rr(:, k+1)=rr(:, k)+(vv(:, k+1)+vv(:, k))*dt/2; %posizione al passo k+1

end

rr_t=rr(:, h); %posizione al tempo t
vv_t=vv(:, h); %velocità al tempo t

figure;
plotObj(6378000, DD, 'b');

if nargin>7
    if Rl=='L'
        plotObj(1737400, dd, 'b');
    else
        plotObj(Rl, dd, 'b');
    end
end
    

rx=rr(1, :);
ry=rr(2, :);
rz=rr(3, :);

plot3(rx, ry, rz, 'k-');

axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
end

