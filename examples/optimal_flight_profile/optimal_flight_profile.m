%% Time Optimal Flight Profile
yops times: t t0 tf
yops state: x size: [7, 1] nominal: [1e3,1e5,1e5,1e2,1,1,1e4]
yops ctrls: alpha phi throttle int: [2,2,0]

flight_param;
d2r    = @(d) d*pi/180;
T_fn   = @(H) Ts + H*dTdH;
p_fn   = @(T) ps*(T/Ts)^(-g/dTdH/R);
rho_fn = @(T) rhos*(T/Ts)^(-(g/dTdH/R+1));
v_fn   = @(H,k) cas2tas( ...
    kappa, p_fn(T_fn(H)), rho_fn(T_fn(H)), ps, rhos, k*ktomps);

H0 = 1600 * ftom;
n0 = 0;
e0 = 0;
v0 = v_fn(H0, 190);
gamma0 = d2r(8);
chi0 = 0;
W0 = 18000 * g;

Hf = 2000*ftom;
nf = 8e5;
ef = 9e5;
vf = v_fn(Hf, 190);
gammaf = d2r(-3);
chif = -3/4*pi;

H_max = 25000 * ftom;
n_max = 12e5;
e_max = 12e5;
v_max = v_fn(H_max, 227);
gamma_max = d2r(60);
chi_max = pi;
W_max = W0;

H_min = H0;
n_min = n0;
e_min = e0;
v_min = v0;
gamma_min = -gamma_max;
chi_min = -pi;
W_min=(18000 - 4000)*g;

x0     = [H0   ; n0   ; e0   ; v0   ; gamma0   ; chi0   ; W0   ];
xf_max = [Hf   ; nf   ; ef   ; vf   ; gammaf   ; chif   ; W0   ];
xf_min = [Hf   ; nf   ; ef   ; vf   ; gammaf   ; chif   ; W_min];
x_max  = [H_max; n_max; e_max; v_max; gamma_max; chi_max; W_max];
x_min  = [H_min; n_min; e_min; v_min; gamma_min; chi_min; W_min]; 

% Control box constraints
u = [alpha; phi; throttle];
du = [der(alpha); der(phi)]; % throttle not rate limited
u_max = [d2r(+10); d2r(+45); 1];
u_min = [d2r(-10); d2r(-45); 0];
du_max = [d2r(+2); d2r(+10)]; % No rate limit on throttle
du_min = [d2r(-2); d2r(-10)];

% Obstacles - circles 
npos = x(2);
epos = x(3);
obs_n = [7.75; 5.45; 5.05; 1.05; 4.45; 6;  1; 8; 3; 6]*1e5; % North
obs_e = [8.42; 4.42; 8.42; 0.72; 0.00; 2;  5; 4; 9; 8]*1e5; % East
obs_r = [   4;    3;    2;    6;    7; 5; 10; 9; 3; 8]*1e4; % radius

%% Solve OCP
ocp = yop.ocp('Optimal Flight Profile');
ocp.min( tf );
ocp.st( t0==0 );
ocp.st( 7e3 <= tf <= 15e3 );
ocp.st( der(x) == flight_dynamics(x, u) );
ocp.st(  x(t0) == x0 );
ocp.st( x_min  <=   x   <= x_max  );
ocp.st( xf_min <= x(tf) <= xf_max );
ocp.st( u_min  <=   u   <= u_max  );
ocp.st( du_min <=  du   <= du_max );

% First solve it without obstacles 
sol = ocp.solve('ival', 100);

% Add the obstacles and resolve the problem with the previous solution as
% initial guess.
ocp.hard( (npos-obs_n).^2 + (epos-obs_e).^2 >= obs_r.^2 );
sol = ocp.solve('ival', 800, 'guess', sol);
%% Plot solution
r2d = @(r) r*180/pi;

figure(1); hold on
for k=1:length(obs_n)
    center=[obs_e(k) obs_n(k)];
    obspos = [center-obs_r(k) 2*obs_r(k) 2*obs_r(k)];
    rectangle('Position',obspos,'Curvature',[1 1]);
end
sol.plot(x(3), x(2))
xlabel('East position [m]'); ylabel('North position [m]')

figure(2); 
subplot(511); hold on
sol.plot(t/3600, x(1))
xlabel('Time [h]'); ylabel('Altitude [m]')
subplot(512); hold on
sol.plot(t/3600, x(4)*3.6)
xlabel('Time [h]'); ylabel('True airspeed [km/h]')
subplot(513); hold on
sol.plot(t/3600, r2d(x(5)))
xlabel('Time [h]'); ylabel('Flight path angle [deg]')
subplot(514); hold on
sol.plot(t/3600, r2d(x(6)))
xlabel('Time [h]'); ylabel('Tracking angle [deg]')
subplot(515); hold on
sol.plot(t/3600, x(7))
xlabel('Time [h]'); ylabel('mass [kg]')

figure(3); 
subplot(311); hold on
sol.plot(t/3600, r2d(u(1)))
xlabel('Time [h]'); ylabel('Angle of attack [deg]')
subplot(312); hold on
sol.plot(t/3600, r2d(u(2)))
xlabel('Time [h]'); ylabel('Roll angle [deg]')
subplot(313); hold on
sol.stairs(t/3600, u(3)*100)
xlabel('Time [h]'); ylabel('Throttle setting [%]')
