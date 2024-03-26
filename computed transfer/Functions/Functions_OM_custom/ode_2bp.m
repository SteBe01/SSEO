function dy = ode_2bp(~, y, mu, R, J, cD, AM, W)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz )     [ L, L/T ]
% mu[1] Gravitational parameter of the primary            [L^3/T^2]
%
% R:  [1] radius of the referenced m1 (usually: Earth)   [km]
% J:  [1] value of the zonal armonic perturbation        [-]
% cD: [1] coefficient of drag                            [-]
% AM: [1] A/M, spaceship's rate of surface over mass     [km^2/kg]
% W : [3x1] planet's rotation angular velocity           [rad/s]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% 2023-11-20: Second version (Pietro Bolsi, Feat: Ancillotti G., Tartaglia D., Tessarollo A.)
%
%
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state

if nargin==3  %unperturbed orbit
    dy = [ v; (-mu/rnorm^3)*r ];
else 
    
    aJ=J*(1.5*mu*R^2/rnorm^5)*[r(1)*((5*r(3)^2/rnorm^2)-1);  %zonal armonic acceleration
                               r(2)*((5*r(3)^2/rnorm^2)-1);
                               r(3)*((5*r(3)^2/rnorm^2)-3)];

    v_rel=v - cross(W, r);

    v_rel_norm = norm(v_rel);

    h=rnorm-R; %[km] actual altitude

    if 0 <= h && h<= 25
        h0 = 0 ;
        rho0 = 1.225*1e9 ; 
        H = 7.249 ; 
    elseif 25 < h && h <= 30
        h0 = 25 ; 
        rho0 = 3.899*1e-2*1e9 ; 
        H = 6.349;
    elseif 30 < h && h <= 40
        h0 = 30 ; 
        rho0 = 1.774*1e-2*1e9 ; 
        H = 6.682;
    elseif 40 < h && h <= 50
        h0 = 40 ; 
        rho0 = 3.972*1e-3*1e9 ; 
        H = 7.554;
    elseif 50 < h && h <= 60
        h0 = 50 ; 
        rho0 = 1.057*1e-3*1e9 ; 
        H = 8.382;
    elseif 60 < h && h <= 70
        h0 = 60 ; 
        rho0 = 3.206*1e-4*1e9 ; 
        H = 7.714;
    elseif 70 < h && h <= 80
        h0 = 70 ; 
        rho0 = 8.770*1e-5*1e9 ; 
        H = 6.549;
    elseif 80 < h && h <= 90
        h0 = 80 ; 
        rho0 = 1.905*1e-5*1e9 ; 
        H = 5.799; 
    elseif 90 < h && h <= 100
        h0 = 90 ; 
        rho0 = 3.396*1e-6*1e9 ; 
        H = 5.382;
    elseif 100 < h && h <= 110
        h0 = 100 ; 
        rho0 = 5.297*1e-7*1e9 ; 
        H = 5.877;
    elseif 110 < h && h <= 120
        h0 = 110 ; 
        rho0 = 9.661*1e-8*1e9 ; 
        H = 7.263;
    elseif 120 < h && h <= 130
        h0 = 120 ; 
        rho0 = 2.438*1e-8*1e9 ; 
        H = 9.473;
    elseif 130 < h && h <= 140
        h0 = 130 ; 
        rho0 = 8.484*1e-9*1e9 ; 
        H = 12.636;
    elseif 140 < h && h <= 150
        h0 = 140 ; 
        rho0 = 3.845*1e-9*1e9 ; 
        H = 16.149;
    elseif 150 < h && h <= 180
        h0 = 150 ; 
        rho0 = 2.070*1e-9*1e9 ; 
        H = 22.523;
    elseif 180 < h && h <= 200
        h0 = 180 ; 
        rho0 = 5.464*1e-10*1e9 ; 
        H = 29.740;
    elseif 200 < h && h <= 250
        h0 = 200 ; 
        rho0 = 2.789*1e-10*1e9 ; 
        H = 37.105;
    elseif 250 < h && h <= 300
        h0 = 250 ; 
        rho0 = 7.248*1e-11*1e9 ; 
        H = 45.546;
    elseif 300 < h && h <= 350
        h0 = 300 ; 
        rho0 = 2.418*1e-11*1e9 ; 
        H = 53.628;
    elseif 350 < h && h <= 400
        h0 = 350 ;
        rho0 = 9.518*1e-12*1e9 ;
        H = 53.298 ; 
    elseif 400 < h && h <= 450
        h0 = 400 ; 
        rho0 = 3.725*1e-12*1e9;
        H = 58.515 ;
    elseif 450 < h && h <= 500 
        h0 = 450 ; 
        rho0 = 1.585*1e-12*1e9 ; 
        H = 60.282 ; 
    elseif 500 < h && h <= 600 
        h0 = 500 ; 
        rho0 = 6.967*1e-13*1e9 ; 
        H = 63.822 ;
    elseif 600 < h && h <= 700 
        h0 = 600 ; 
        rho0 = 1.454*1e-13*1e9 ; 
        H = 71.835 ;
    elseif 700 < h && h <= 800 
        h0 = 700 ; 
        rho0 = 3.614*1e-14*1e9 ; 
        H = 88.667 ;
    elseif 800 < h && h <= 900 
        h0 = 800 ; 
        rho0 = 1.170*1e-14*1e9 ; 
        H = 124.64 ; 
    elseif 900 < h && h <= 1000 
        h0 = 900 ; 
        rho0 =5.245*1e-15*1e9 ; 
        H = 181.05 ; 
    elseif h > 1000  
        h0 = 1000 ; 
        rho0 = 3.019*1e-15*1e9 ; 
        H = 268 ;
    end

    rho=rho0*exp((h0-h)/H);

    aAD=-0.5*rho*AM*cD*v_rel_norm*(v_rel);  %air drag acceleration

    dy = [ v; ((-mu/rnorm^3)*r) + aJ + aAD];

end
end
