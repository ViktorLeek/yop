yops Times: t t0 tf
yops States: rad lon lat v fpa azi nominal: [1e6,1,1,1e4,1,1]
yops Controls: aoa bank

x = [rad; lon; lat; v; fpa; azi]; 
u = [aoa; bank]; 

d2r   = pi/180;
alt0  = 79248;
Re    = 6371204; % Earth radius [m]
x0    = [alt0+Re;   0;       0;  7803;  -1*d2r;   90*d2r];
x_max = [alt0+Re;  pi;  70*d2r; 45000;  80*d2r;  180*d2r];
x_min = [     Re; -pi; -70*d2r;    10; -80*d2r; -180*d2r];
u_max = [ 90*d2r;   1*d2r];
u_min = [-90*d2r; -90*d2r];

%% First guess
sim = yop.ivp(t0==0, tf==2000);
sim.add( der(x) == launch_vehicle(x, u) );
sim.add(  x(t0) == x0 );
sim.add(  u(t)  == 0 );
res = sim.solve('points', 2000);

%% Plot guess
figure(1)
subplot(321); hold on;
res.plot(t, rad - Re)
subplot(322); hold on;
res.plot(t, v)
subplot(323); hold on;
res.plot(lon*180/pi, lat*180/pi)
subplot(324); hold on;
res.plot(t, fpa*180/pi)
subplot(325); hold on;
res.plot(t, aoa*180/pi)
subplot(326); hold on;
res.plot(t, bank*180/pi)

%% Reusable Launch Vehicle Entry
ocp = yop.ocp('Reusable Launch Vehicle Entry');
ocp.max( lat(tf) );
ocp.st( t0 == 0 );
ocp.st( 1000 <= tf <= 3000 );
ocp.st( der(x) == launch_vehicle(x, u) );
ocp.st(  x(t0) == x0                   );
ocp.st( x_min <= x <= x_max );
ocp.st( u_min <= u <= u_max );
ocp.st( rad(tf) == 24384 + Re );
ocp.st(   v(tf) == 762        );
ocp.st( fpa(tf) == -5*d2r     );

sol = ocp.solve('guess', res, 'ival', 250, 'dx', 3);

%% Plot solution
figure(1)
subplot(321); hold on;
sol.plot(t, (rad - Re)*1e-3)
xlabel('t [s]')
ylabel('h [km]')

subplot(322); hold on;
sol.plot(t, v*1e-3)
xlabel('t [s]')
ylabel('v [km/s]')

subplot(323); hold on;
sol.plot(lon*180/pi, lat*180/pi)
xlabel('\phi [deg]')
ylabel('\theta [deg]')

subplot(324); hold on;
sol.plot(t, fpa*180/pi)
xlabel('t [s]')
ylabel('\gamma [deg]')

subplot(325); hold on;
sol.plot(t, aoa*180/pi)
xlabel('t [s]')
ylabel('\alpha [deg]')

subplot(326); hold on;
sol.plot(t, bank*180/pi)
xlabel('t [s]')
ylabel('\sigma [deg]')

%% Model
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
cl = [-0.2070; 1.6756];         % Lift coeffcicient
cd = [0.0785; -0.3529; 2.0400]; % Drag coefficient
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
