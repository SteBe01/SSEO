function dx = ode_3bp(~,x,MU,mu)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dx = ode_3bp(t,x,MU,mu)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% x[12x1] State of the Earth and of the body ( rx, ry, rz, vx, vy, vz, rx, ry, rz, vx, vy, vz )     [ L, L/T, L, L/T ] from the Sun
% MU[1] Gravitational parameter of the Sun             [L^3/T^2]
% mu[1] Gravitational parameter of the Earth           [L^3/T^2]
%
% OUTPUT:
% dx[12x1] Derivative of the state [ L/T^2, L/T^3, L/T^2, L/T^3  ]
%
% CONTRIBUTORS:
% Pietro Bolsi

% Position and velocity of the Earth from the Sun
R = x(1:3);
V = x(4:6);
% Position and velocity of the body from the Sun
r = x(7:9);
v = x(10:12);

% Distance of Earth from Sun
n1 = norm(R);
% Distance of body from Sun
n2 = norm(r);
% Distance of body from Earth
n3 = norm(R-r);

% Set the derivatives of the state

dx = [ V; (-MU/n1^3)*R; v; (-MU/n2^3)*r + (-mu/n3^3)*(r-R) ];

end


