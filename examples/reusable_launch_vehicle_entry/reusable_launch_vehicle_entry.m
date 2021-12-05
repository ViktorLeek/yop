clear;

t  = yop.time();
t0 = yop.time0();
tf = yop.timef();

rad  = yop.state();
lon  = yop.state();
lat  = yop.state();
v    = yop.state();
fpa  = yop.state();
azi  = yop.state();
x = [rad; lon; lat; v; fpa; azi]; 

aoa  = yop.control;%('pw', 'quadratic');
bank = yop.control;%('pw', 'linear');
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
ivp = yop.ivp( ...
    t0==0, tf==2000, ...
    der(x) == launch_vehicle(x, u), ...
    x(t0)  == x0, ...
    u == [17.4*d2r; 70/2000*t-70] ...
    );
sim = ivp.solve();

%% Improved guess (relaxed problem - fixed time)
ocp = yop.ocp('Reusable Launch Vehicle Entry - fixed time');
ocp.max(lat(tf));
ocp.st(...
    tf == 2000, ...
    der(x) == launch_vehicle(x, u), ...
     x(t0) == x0, ...
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    rad(tf) == 24384 + Re, ...
      v(tf) == 762, ...
    fpa(tf) == -5*d2r, ...
    azi(tf) == -90*d2r ...
    );

sol = ocp.solve('intervals', 200, 'degree', 3, 'guess', sim);
% sol.plot(t, rad - Re)
%% Reusable Launch Vehicle Entry
ocp = yop.ocp('Reusable Launch Vehicle Entry');
ocp.max(lat(tf));
ocp.st(...
    1000 <= tf <= 3000, ...
    der(x) == launch_vehicle(x, u), ...
     x(t0) == x0, ...
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    rad(tf) == 24384 + Re, ...
      v(tf) == 762, ...
    fpa(tf) == -5*d2r, ...
    azi(tf) == -90*d2r ...
    );

sol2 = ocp.solve('intervals', 1000, 'degree', 5, 'guess', sol);


%%
figure(1)
subplot(321); hold on;
sol2.plot(t, rad - Re)
subplot(322); hold on;
sol2.plot(t, v)
subplot(323); hold on;
sol2.plot(lon*180/pi, lat*180/pi)
subplot(324); hold on;
sol2.plot(t, fpa*180/pi)
subplot(325); hold on;
sol2.plot(t, aoa*180/pi)
subplot(326); hold on;
sol2.plot(t, bank*180/pi)
















