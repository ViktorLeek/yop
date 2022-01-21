%% Second order singular regulator
%   From: https://www.bocop.org/regulator/
yops Times: t t0 tf 
yops States: x size: [2,1]
yops Control: u

ocp = yop.ocp('Fuller Problem');
ocp.min( int( x(1)^2 + der(x(1))^2 ) );
ocp.st( t0==0, tf==5 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [   0; 1] );
ocp.st( -1 <= u <= 1  );
sol = ocp.solve('ival', 1000, 'dx', 2);

figure(1); 
subplot(311); hold on;
sol.plot(t, x(1))
subplot(312); hold on;
sol.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)
