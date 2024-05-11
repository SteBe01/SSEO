function [dv, rp, v1p, v2p] = fly_by(vI, vF, IDplanet)
%INPUTs:
% vI       : incoming velocity vector    [3X1]    [km/s]
% vF       : outgoing velocity vector    [3X1]    [km/s]
% IDplanet : planet ID                   [1X1]    [-]
%
%OUTPUTs:
% dv   : delta velocity for the flyby at pericentre   [1X1]  [km/s]
% rp   : pericentre radius                            [1X1]  [km]
% v1p  : ingoing velocity at pericentre               [1X1]  [km/s]
% v2p  : outgoing velocity at pericentre              [1X1]  [km/s]
%
% if rp is less than the radius of the planet +150 km, the function gives
% NaN as output for dv, v1p, v2p.
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

mu=astroConstants(10+ IDplanet); 
rp_lim=astroConstants(20+ IDplanet) + 150; %[km] radius of the planet
rp0=[0, 1e5*rp_lim];%[km]


if size(vI)==[1, 3]  %the function works with column vectors only!
    vI=vI';
end

if size(vF)==[1, 3]  %the function works with column vectors only!
    vF=vF';
end

v1=norm(vI);       %incolming velocity (modulus)
v2=norm(vF);       %exit velocity (modulus)

format long

d=real(acos(dot(vI, vF)/(v1*v2)));    %delta

if d<=0
    rp=Inf;

    v1p=v1;
    v2p=v2;

    dv=abs(v2-v1);
elseif d>=pi
    rp=0;

    v1p=NaN;
    v2p=NaN;

    dv=NaN;
else
    if 1==isnan(v1)
        rp=NaN;
        v1p=NaN;
        v2p=NaN;
        dv=NaN;
    elseif 1==isnan(v2)
        rp=NaN;
        v1p=NaN;
        v2p=NaN;
        dv=NaN;
    else
        e1=@(x) 1+x.*(v1^2)/mu;       %incolming eccentricity

        e2=@(y) 1+y.*(v2^2)/mu;       %exit eccentricity
    
        f=@(z) real(asin(1/e1(z)))+real(asin(1/e2(z)))-d;

        rp=fzero(f, rp0);    %find the radius at pericentre
    
        v1p=sqrt(v1^2 + 2*mu/rp);  %incoming velocity at pericentre
        v2p=sqrt(v2^2 + 2*mu/rp);  %exit velocity at pericentre

        dv=abs(v2p-v1p);  %delta velocity for the flyby at pericentre
    end
end
 
if astroConstants(13)*0.99<mu && mu<astroConstants(13)*1.01 %astroConstants(13) %mu=mu_earth
    if rp<=rp_lim
        %warning("The satellite is colliding with the Earth!");
        %rp=NaN;
        v1p=NaN;
        v2p=NaN;
        dv=NaN;
    elseif 1==isnan(rp)
        disp(d);
        v1p=NaN;
        v2p=NaN;
        dv=NaN;
    end


end

end

