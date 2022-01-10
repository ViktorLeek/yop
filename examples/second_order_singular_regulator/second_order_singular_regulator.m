%% Second order singular regulator
%   From: https://www.bocop.org/regulator/
yops times: t t0 tf 
yops states: x size: [2,1]
yops control: u

ocp = yop.ocp('Fuller Problem');
ocp.min( int( x(1)^2 + der(x(1))^2 ) );
ocp.st( tf == 5 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [   0; 1] );
ocp.st( -1 <= u <= 1  );
sol = ocp.solve('intervals', 1000, 'degree', 2);

figure(1); 
subplot(311); hold on;
sol.plot(t, x(1))
subplot(312); hold on;
sol.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)
