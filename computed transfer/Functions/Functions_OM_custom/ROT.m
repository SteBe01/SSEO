function R = ROT(i,OM,om)
% this function creates the rotation matrix to go from the perifocal to
% eci/heliocentric frame
% 
%INPUTS:
% i  : inclination of the orbit        [1x1] [rad]
% OM : RAAN of the orbit               [1x1] [rad]
% om : pericentre anomaly of the orbit [1x1] [rad]
%
%OUTPUT:
% R  : rotation matrix for the position and velocity vectors [3x3] [-]
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

    ROM = [cos(OM) sin(OM)   0  ;-sin(OM) cos(OM)     0  ;  0      0       1  ];
    Ri  = [  1        0      0  ;    0    cos(i)   sin(i);  0   -sin(i) cos(i)];
    Rom = [cos(om) sin(om)   0  ;-sin(om) cos(om)     0  ;  0      0       1  ];
    R=ROM'*Ri'*Rom';
end


