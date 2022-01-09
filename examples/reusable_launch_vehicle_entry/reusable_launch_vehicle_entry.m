yops times: t t0 tf
yops states: rad lon lat v fpa azi scaling: [1e6,1,1,1e4,1,1]
yops ctrls: aoa bank

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
ocp.st( 1000 <= tf <= 3000 );
ocp.st( der(x) == launch_vehicle(x, u) );
ocp.st(  x(t0) == x0                   );
ocp.st( x_min <= x <= x_max );
ocp.st( u_min <= u <= u_max );
ocp.st( rad(tf) == 24384 + Re );
ocp.st(   v(tf) == 762        );
ocp.st( fpa(tf) == -5*d2r     );

sol = ocp.solve('guess', res, 'intervals', 250, 'degree', 3);

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
