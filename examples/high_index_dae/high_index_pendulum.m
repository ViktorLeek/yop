%% Variables and parameters
% Problem from: 
%   Campbell, Kunkel - Solving higher index DAE optimal control problems
%   https://www.aimsciences.org/article/doi/10.3934/naco.2016020
yops Times: t t0 tf
yops States: x1 x2 x3 x4
yops Control: u
yops Algebraic: F 
yops Control: slack

x = [x1; x2; x3; x4];
a = 0.5; 
c = 1;
d = 100;
g = 4;
L = 2;
alpha = 1*pi/180;

%% OCP - Index 3, No index reduction or stabilization (Fails, with mumps)
ocp = yop.ocp();
ocp.min( int( c*u^2 + d*(x1 - L*sin(t+alpha))^2 + d*(x3 - L*cos(t+alpha))^2 ));
ocp.st( t0==0, tf==2.2 );
ocp.st( der(x1) == x2 );
ocp.st( der(x2) == -F*x1 - a*x2 + u*x3 );
ocp.st( der(x3) == x4 );
ocp.st( der(x4) == -F*x3 - a*x4 - g -u*x1 );
ocp.st( x(t0) == [0; 0; L; 0] );
ocp.alg( 0 == x1^2 + x3^2 - L^2);
sol = ocp.solve();

%% OCP - Index 3, relaxed problem as initial guess (converges to solution)
rxd = yop.ocp();
rxd.min( int( c*u^2 + d*(x1 - L*sin(t+alpha))^2 + d*(x3 - L*cos(t+alpha))^2 ) + 1e5*int(slack^2));
rxd.st( t0==0, tf==2.2 );
rxd.st( der(x1) == x2 );
rxd.st( der(x2) == -F*x1 - a*x2 + u*x3 );
rxd.st( der(x3) == x4 );
rxd.st( der(x4) == -F*x3 - a*x4 - g -u*x1 );
rxd.st( x(t0) == [0; 0; L; 0] );
rxd.alg( 0 == x1^2 + x3^2 - L^2 + slack);
rsl = rxd.solve();
sol = ocp.solve('guess', rsl);

figure(3); rsl.plot(t, slack);

%% OCP - Overdetermined (Not solvable, with IPOPT)
ocp = yop.ocp();
ocp.min( int( c*u^2 + d*(x1 - L*sin(t+alpha))^2 + d*(x3 - L*cos(t+alpha))^2 ) );
ocp.st( t0==0, tf==2.2 );
ocp.st( der(x1) == x2 );
ocp.st( der(x2) == -F*x1 - a*x2 + u*x3 );
ocp.st( der(x3) == x4 );
ocp.st( der(x4) == -F*x3 - a*x4 - g -u*x1 );
ocp.st( x(t0) == [0; 0; L; 0] );
ocp.alg( 0 == x1^2 + x3^2 - L^2);
ocp.alg( 0 == x1*x2 + x3*x4 );
ocp.alg( 0 == x2^2 + x4^2 -F*L^2 - x3*g );
sol = ocp.solve();

%% Semi-explicit index-1 (More solutions than original problem)
ocp = yop.ocp();
ocp.min( int( c*u^2 + d*(x1 - L*sin(t+alpha))^2 + d*(x3 - L*cos(t+alpha))^2 ) );
ocp.st( t0==0, tf==2.2 );
ocp.st( der(x1) == x2 );
ocp.st( der(x2) == -F*x1 - a*x2 + u*x3 );
ocp.st( der(x3) == x4 );
ocp.st( der(x4) == -F*x3 - a*x4 - g -u*x1 );
ocp.st( x(t0) == [0; 0; L; 0] );
ocp.alg( 0 == x2^2 + x4^2 -F*L^2 - x3*g );
sol = ocp.solve();

%% Baumgarte stabilized OCP
k = -15;
b = 50;
ocp = yop.ocp();
ocp.min( int( c*u^2 + d*(x1 - L*sin(t+alpha))^2 + d*(x3 - L*cos(t+alpha))^2 ) );
ocp.st( t0==0, tf==2.2 );
ocp.st( der(x1) == x2 );
ocp.st( der(x2) == -F*x1 - a*x2 + u*x3 );
ocp.st( der(x3) == x4 );
ocp.st( der(x4) == -F*x3 - a*x4 - g -u*x1 );
ocp.st( x(t0) == [0; 0; L; 0] );
ocp.alg( 0 == (x2^2 + x4^2 -F*L^2 - x3*g) + k*(x1*x2 + x3*x4) + b*(x1^2 + x3^2 - L^2) );
sol = ocp.solve();

%% Residuals
figure(1); 
subplot(311); hold on
sol.plot(t, x1^2 + x3^2 - L^2)
subplot(312); hold on
sol.plot(t, x1*x2 + x3*x4)
subplot(313); hold on
sol.plot(t, x2^2 + x4^2 -F*L^2 - x3*g)

%% Pendulum path and control
figure(2); 
subplot(211); hold on
sol.plot(x1, x3)
subplot(212); hold on
sol.stairs(t, u)