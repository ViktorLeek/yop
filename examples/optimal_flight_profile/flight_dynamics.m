function dx = flight_dynamics(x,u)
% Aircraft Dynamics - Internal
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%    g_eq - constraint function for equality constraints
%    g_neq - constraint function for inequality constraints
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

H     = x(1);
npos  = x(2);
epos  = x(3);
v     = x(4); % TAS
gamma = x(5);
chi   = x(6);
W     = x(7);
alpha    = u(1);
phi      = u(2);
throttle = u(3);

flight_param;

% Look up tables
a1p =  1.5456e-9;
a2p = -3.1176e-7;
a3p =  1.9477e-5;
b1p = -2.0930e-5;
b2p =  4.1670e-3;
b3p = -4.0739e-1;
c1p =  8.8835e-2;
c2p = -1.2855e+1;
c3p =  3.0777e+3;
Pidle = 75;

Temp = Ts + H*dTdH;
pressure = ps * (Temp/Ts) ^ (-g / dTdH / R);
rho = rhos * (Temp/Ts) ^ (-(g / dTdH / R + 1));
V_cas = cas2tas(kappa, pressure, rho, ps, rhos, v);

ap = a1p*V_cas^2 + a2p*V_cas + a3p;
bp = b1p*V_cas^2 + b2p*V_cas + b3p;
cp = c1p*V_cas^2 + c2p*V_cas + c3p;
Pmax = ap*H^2 + bp*H + cp;

P = (Pmax - Pidle)*throttle + Pidle;
J = 60 * v / nprop / Dprop;
ita = -0.13289*J^6 + 1.2536*J^5 - 4.8906*J^4 + 10.146*J^3 - 11.918*J^2 ...
    + 7.6740*J - 1.3452;
Thrust = 2*745.6*ita*P/v;

% Aerodynamic forces
cl   = clalpha * (alpha - alpha0) * 180/pi;
cd   = cd0 + k_cd * cl^2;
L    = 0.5 * rho * S * cl * v^2;
Drag = 0.5 * rho * S * cd * v^2;

% Equations of motions
v_dot = ( Thrust - Drag - W*sin(gamma) ) * g / W;
gamma_dot = ( L*cos(phi) + Thrust*sin(alphat) - W*cos(gamma) )*g/W/v;
H_dot = v * sin(gamma);
x_dot = v * cos(gamma) * cos(chi);
y_dot = v * cos(gamma) * sin(chi);
chi_dot = L * sin(phi) / cos(gamma) * g / W / v;
W_dot = weight_rate(V_cas, H, throttle);
dx = [H_dot; x_dot; y_dot; v_dot; gamma_dot; chi_dot; W_dot];
end

function dw = weight_rate(v_cas, H, throttle)

% coefficients for maximum throttle setting
am1 = -1.5374e-13;
am2 = +2.3875e-11;
am3 = -5.7989e-10;
bm1 = +2.1483e-09;
bm2 = -3.2731e-07;
bm3 = -9.1353e-07;
cm1 = -5.8602e-06;
cm2 = +9.7605e-04;
cm3 = +1.1389e-01; 

% coefficients for minimum throttle setting (i for idle)
ai1 = +4.2242e-15;        
ai2 = -9.0868e-13;       
ai3 = +1.6801e-10;    
bi1 = -4.8122e-11;   
bi2 = +1.0818e-08;     
bi3 = -2.7199e-06;    
ci1 = +2.5788e-08;     
ci2 = -6.8894e-06;    
ci3 = +2.5804e-02;        
  
% Maximum fuel flow
am = am1*v_cas^2 + am2*v_cas + am3;
bm = bm1*v_cas^2 + bm2*v_cas + bm3;
cm = cm1*v_cas^2 + cm2*v_cas + cm3;
Fmax = am*H^2 + bm*H + cm; % maximum fuel flow    [kg/s]

% Minimum fuel flow (idle throttle setting)
ai = ai1*v_cas^2 + ai2*v_cas + ai3;
bi = bi1*v_cas^2 + bi2*v_cas + bi3;
ci = ci1*v_cas^2 + ci2*v_cas + ci3;
F_idle = ai*H^2 + bi*H + ci; % minimum fuel flow [kg/s]

% Interpolation
u = 0.08*throttle.^2 + 0.92*throttle; % Actuator map
FF = F_idle + (Fmax - F_idle)*u; % Current fuel flow [kg/s]

% Calculation of weight
dw = -2*FF*9.81;
end


