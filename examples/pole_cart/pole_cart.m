%yops Times: t t0 tf State: x size: [4,1] Control: u
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state('size', [4,1]);
u  = yop.control();

p = x(1); % Cart position
v = x(2); % Cart speed
o = x(3); % Pole angle
w = x(4); % Pole angular velocity

Q = diag([1e3, 1e-2, 1e3, 1e-2]);
R = 1e-2;
J = x'*Q*x + u'*R*u;
ocp = yop.ocp();
ocp.min( J(tf) + int(J) );
ocp.st( t0==0, tf==3 );
ocp.st( der(x) == cart(x,u) );
ocp.st(  x(t0) == [0; 0; pi; 0]);
ocp.st( -15 <= u <= 15 );

sol = ocp.solve('ival', 100, 'dx', 3);

figure(1); hold on
sol.stairs(t, u/10)
sol.plot(t, x(1:3))

function dx = cart(x,u)
M = 1;    % Cart mass [kg]
m = 0.1;  % Ball mass [kg]
l = 0.8;  % Rod length [m]
g = 9.81; % gravitional acceleration [m/s^2]

% u = Cart force
p = x(1); % Cart position
v = x(2); % Cart speed
o = x(3); % Pole angle
w = x(4); % Pole angular velocity

dv = (-l*m*sin(o)*w^2 + u + g*m*cos(o)*sin(o))/(M + m - m*cos(o)^2);
dw = (-l*m*cos(o)*sin(o)*w^2 + u*cos(o) + (m+M)*g*sin(o))/(l*M + l*m*(1-cos(o)^2));
dx = [v;dv;w;dw];
end