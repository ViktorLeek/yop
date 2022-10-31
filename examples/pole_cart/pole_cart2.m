t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state('size', [4,1]);
u  = yop.control('int', 1);

Q = diag([10, 2.5, 1, 0.1]);

ocp = yop.ocp();
ocp.min( int(u^2 + x'*Q*x) );
ocp.st( t0==0, tf==2 );
ocp.st( der(x) == polecart(x,u) );
ocp.st( x(t0) == 0 );
ocp.st( x(tf) == [1; 0; pi; 0] );
ocp.st( -2 <= x(1) <=  2 );
ocp.st( -20 <=  u   <= 20 );
ocp.hard( 20 <= der(0+x(4)) <= 20 );

sol = ocp.solve('ival', 100, 'dx', 3);

figure(1); hold on
sol.plot(t, x)
sol.plot(t, u)

function dx = polecart(x,u)
M = 1;    % Cart mass [kg]
m = 0.3;  % Ball mass [kg]
l = 0.5;  % Pole length [m]
g = 9.81; % gravitional acceleration [m/s^2]

p = x(1); % Cart position
v = x(2); % Cart speed
o = x(3); % Pole angle
w = x(4); % Pole angular velocity

dv =  (l*m*sin(o)*w^2+u+m*g*cos(o)*sin(o))/(M+m*(1-cos(o)^2));
dw = -(l*m*cos(o)*sin(o)*w^2+u*cos(o)+(m+M)*g*sin(o))/(l*M+l*m*(1-cos(o)^2));
dx = [v;dv;w;dw];
end