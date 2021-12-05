function dx = launch_vehicle(x, u)

rad = x(1); 
lon = x(2); 
lat = x(3); 
v   = x(4); 
fpa = x(5); 
azi = x(6);

aoa  = u(1); 
bank = u(2); 

Re = 6371204;                   % Earth radius [m]
S  = 249.1;                     % Vehicle reference area [m2]
cl = [-0.2070, 1.6756];         % Lift coeffcicient
cd = [0.0785, -0.3529, 2.0400]; % Drag coefficient
H  = 7254.24;                   % Density scale height [m]
rho0 = 1.225570827014494;       % Density of air at sea level [kg/m3]
mu   = 3.986031954093051e14;    % Earth gravitational paramter [m3/s2]
mass = 92079.2525560557;        % Vehicle mass [kg]

alt = rad - Re;
cD  = cd(1) + cd(2)*aoa + cd(3)*aoa^2;
rho = rho0*exp(-alt/H);
cL  = cl(1) + cl(2)*aoa;
q   = 0.5*rho*v^2;
D   = q*S*cD/mass;
L   = q*S*cL/mass;
g   = mu/rad^2;

dx = [ ...
    v*sin(fpa); ...
    v*cos(fpa)*sin(azi)/(rad*cos(lat)); ...
    v*cos(fpa)*cos(azi)/rad; ...
    -D - g*sin(fpa); ...
    ( L*cos(bank) - cos(fpa)*(g-v^2/rad) )/v; ...
    ( L*sin(bank)/cos(fpa) + v^2*cos(fpa)*sin(azi)*tan(lat)/rad )/v ...
    ];
end