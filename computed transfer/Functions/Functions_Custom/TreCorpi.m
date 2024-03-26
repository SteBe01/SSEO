function [rr_t, vv_t] = TreCorpi(M, m, d, rr0, vv0, t, dt, Rl)
%TERRALUNA Summary of this function goes here
%   Detailed explanation goes here
G=6.6743e-11; %[m^3/(kg*s^2)]

if nargin<6

    dt=10;

end

X=[1; 0; 0];
%Y=[0; 1; 0];
%Z=[0; 0; 1];

DD=-X*(d*m/(M+m)); %posizione di M rispetto a O all'istante t=0s
dd=+X*(d*M/(M+m)); %posizione di m rispetto a O all'istante t=0s

om=sqrt((G*(M+m))/(d^3));%velocità angolare angolare del sistema omm=[0, 0, om];
omm=[0, 0, om];

tt=[0:dt:t];%vettore dei tempi
h=length(tt);



rr=zeros(3, h); %inizializzo posizioni del satellite
vv=zeros(3, h); %inizializzo velocità (plurale) satellite
aag=zeros(3, h); %inizializzo accelerazioni satellite
dr=zeros(3, h); %vettorini traiettoria del satellite
dv=zeros(3, h);

ppm=zeros(3, h);
ppM=zeros(3, h);

rrm=zeros(3, h);
rrM=zeros(3, h);
%rt=zeros(3, h);
aam=zeros(3, h);
aaM=zeros(3, h);
%vr=zeros(3, h);
rrs=zeros(3, h);

vvt=cross(omm, DD)'; %velocità di M rispetto a O al tempo t=0

rr(:, 1)=rr0+DD; %posizione satellite rispetto a O al tempo t=0 (riferito al sistema fisso)
vv(:, 1)=vv0+vvt; %velocità del satellite rispetto a O al tempo t=0 (riferita al sistema fisso)



%vr(:, 1)=vv0; %velocità relativa del satellite rispetto a M

for k=1:h

    ROT=[cos(om*tt(k)), sin(om*tt(k)), 0; -sin(om*tt(k)), cos(om*tt(k)), 0; 0, 0, 1]; %matrice di rotazione di del sistema al passo k
    ppm(:, k)=(ROT')*dd;     %posizione di m al passo k 
    ppM(:, k)=(ROT')*DD;     %posizione di M al passo k 

    rrm(:, k)=rr(:, k)-ppm(:, k); %posizione al passo k del satellite rispetto a m
    rrM(:, k)=rr(:, k)-ppM(:, k); %posizione al passo k del satellite rispetto a M

    aam(:, k)=-G*m*rrm(:, k)/(norm(rrm(:, k))^3); %accelerazione grav per m al passo k
    aaM(:, k)=-G*M*rrM(:, k)/(norm(rrM(:, k))^3); %accelerazione grav per M al passo k

    aag(:, k)=aam(:, k)+aaM(:, k); %accelerazione complessiva al passo k

    dv(:, k)=aag(:, k)*dt; %delta di velocità tra passo k+1 e k

    vv(:, k+1)=vv(:, k)+dv(:, k); %velocità al passo k+1
    %vr(:, k+1)=vr(:, k)+dv(:, k); %velocità relativa al psso k+1

    dr(:, k)=(vv(:, k+1)+vv(:, k))*dt/2; %delta spostamento al passo k

    rr(:, k+1)=rr(:, k)+dr(:, k); %posizione al passo k+1

    rrs(:, k)=ROT*rr(:, k);

end

% plotto sistema stelle fisse

rx=rr(1, :);
ry=rr(2, :);
rz=rr(3, :);

pMx=ppM(1, :);
pMy=ppM(2, :);
pMz=ppM(3, :);

figure;

plot3(rx, ry, rz, 'm-', pMx, pMy, pMz, 'y.');
hold on;

if nargin>7

    pmx=ppm(1, :);
    pmy=ppm(2, :);
    pmz=ppm(3, :);

    plot3(pmx, pmy, pmz, 'g.');

end

axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

% plotto sistema rotante

rsx=rrs(1, :);
rsy=rrs(2, :);
rsz=rrs(3, :);

figure;

plotObj(6378000, DD, 'b');

if nargin>7
    if Rl=='L'
        plotObj(1737400, dd, 'b');
    else
        plotObj(Rl, dd, 'b');
    end
end
    
plot3(rsx, rsy, rsz, 'k-');

axis equal;
xlabel('X"');
ylabel('Y"');
zlabel('Z');

rr_t=rrs(:, h); %posizione al tempo t
%vv_t=vvs(:, h); %velocità al tempo t
end

