function [dx, v_vec, Qr, r_vec] = spacecraft(x, u, tau)

p=x(1); f=x(2); g=x(3); h=x(4); k=x(5); L=x(6); w=x(7);
ur=u(1); ut=u(2); uh=u(3);

T   = 4.446618e-3;    % [lb]
Isp = 450;            % [s]
mu  = 1.407645794e16; % [ft^3/s^2]
gs  = 32.174;         % [ft/s^2]
Re  = 20925662.73;    % [ft]
J2  = 1082.639e-6;
J3  = -2.565e-6;
J4  = -1.608e-6;

% Gravitational disturbing acceleration
q   = 1 + f*cos(L) + g*sin(L);
r   = p/q;
a2  = h^2 - k^2;
chi = sqrt(h^2 + k^2);
s2  = 1 + chi^2;

rx = (r/s2) * (cos(L) + a2*cos(L) + 2*h*k*sin(L));
ry = (r/s2) * (sin(L) - a2*sin(L) + 2*h*k*cos(L));
rz = (2*r/s2) * (h*sin(L) - k*cos(L));
r_vec = [rx; ry; rz];
r_mag = norm(r_vec);

vx = -(1/s2) * sqrt(mu/p) * (  sin(L) + a2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + a2*g);
vy = -(1/s2) * sqrt(mu/p) * ( -cos(L) + a2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + a2*f);
vz =  (2/s2) * sqrt(mu/p) * (h*cos(L) +  k*sin(L) + f*h + g*k);
v_vec = [vx; vy; vz];

r_Xv = cross(r_vec, v_vec);
r_Xv_mag = norm(r_Xv);
r_Xv_Xr  = cross(r_Xv, r_vec);

ir1 = rx/r_mag;
ir2 = ry/r_mag;
ir3 = rz/r_mag;
ir = [ir1; ir2; ir3];
it1 = r_Xv_Xr(1) / (r_Xv_mag * r_mag);
it2 = r_Xv_Xr(2) / (r_Xv_mag * r_mag);
it3 = r_Xv_Xr(3) / (r_Xv_mag * r_mag);
it = [it1; it2; it3];
ih1 = r_Xv(1) / r_Xv_mag;
ih2 = r_Xv(2) / r_Xv_mag;
ih3 = r_Xv(3) / r_Xv_mag;
ih = [ih1; ih2; ih3];
enir = ir3;

enirir1 = enir * ir1;
enirir2 = enir * ir2;
enirir3 = enir * ir3;
enenirir1 = 0 - enirir1;
enenirir2 = 0 - enirir2;
enenirir3 = 1 - enirir3;
enenirir_mag = norm([enenirir1, enenirir2, enenirir3]);
in1 = enenirir1 / enenirir_mag;
in2 = enenirir2 / enenirir_mag;
in3 = enenirir3 / enenirir_mag;

% Geocentric latitude
sinphi = rz/norm([rx, rz]);
cosphi = sqrt(1 - sinphi^2);

% Legendre polynomials
P2  = 0.5000 * ( 3 * sinphi^2 - 2);
P3  = 0.5000 * ( 5 * sinphi^3 - 3*sinphi);
P4  = 0.1250 * (35 * sinphi^4 - 30*sinphi^2 + 3);
dP2 = 3.0000 * sinphi;
dP3 = 0.5000 * ( 15 * sinphi - 3);
dP4 = 0.1250 * (140 * sinphi^3 - 60*sinphi);


% Oblate earth perturbations
sumn =   (Re/r)^2 * dP2*J2 +   (Re/r)^3 * dP3*J3 +   (Re/r)^4 * dP4*J4;
sumr = 3*(Re/r)^2 *  P2*J2 + 4*(Re/r)^3 *  P3*J3 + 5*(Re/r)^4 *  P4*J4; 
Dgn = -mu*cosphi/(r^2) * sumn;
Dgr = -mu/(r^2)*sumr;
Dgnin1 = Dgn * in1; 
Dgnin2 = Dgn * in2; 
Dgnin3 = Dgn * in3;
Dgrir1 = Dgr * ir1; 
Dgrir2 = Dgr * ir2; 
Dgrir3 = Dgr * ir3;
Dg1 = Dgnin1 - Dgrir1; 
Dg2 = Dgnin2 - Dgrir2; 
Dg3 = Dgnin3 - Dgrir3;
Deltag1 = ir(1)*Dg1 + ir(2)*Dg2 + ir(3)*Dg3; 
Deltag2 = it(1)*Dg1 + it(2)*Dg2 + it(3)*Dg3; 
Deltag3 = ih(1)*Dg1 + ih(2)*Dg2 + ih(3)*Dg3;

% Thrust acceleration
DeltaT1 = gs*T*(1 + 0.01*tau)/w * ur;
DeltaT2 = gs*T*(1 + 0.01*tau)/w * ut;
DeltaT3 = gs*T*(1 + 0.01*tau)/w * uh;

% Total acceleration
Delta1 = Deltag1 + DeltaT1;
Delta2 = Deltag2 + DeltaT2;
Delta3 = Deltag3 + DeltaT3;

% Differential equations of motion
dp = 2*p/q * sqrt(p/mu) * Delta2;
df =  sqrt(p/mu) * sin(L)*Delta1 + sqrt(p/mu)*(1/q)*((q+1)*cos(L) + f)*Delta2 - sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L))*Delta3;
dg = -sqrt(p/mu) * cos(L)*Delta1 + sqrt(p/mu)*(1/q)*((q+1)*sin(L) + g)*Delta2 + sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L))*Delta3;
dh =  sqrt(p/mu) * (s2*cos(L)/(2*q))*Delta3;
dk =  sqrt(p/mu) * (s2*sin(L)/(2*q))*Delta3;
dL =  sqrt(p/mu) * (1/q)*(h*sin(L) - k*cos(L))*Delta3 + sqrt(mu*p)*((q/p)^2);
dw = -(T*(1+0.01*tau)/Isp);
dx = [dp; df; dg; dh; dk; dL; dw];
Qr = [ir, it, ih];

