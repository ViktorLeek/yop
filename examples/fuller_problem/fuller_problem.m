%% Fuller Problem
%   From: https://www.bocop.org/fuller-problem/
yops Times: t t0 tf 
yops States: x size: [2,1] scaling: [10,1]
yops Control: u scaling: 1e-2

ocp = yop.ocp('Fuller Problem');
ocp.min( int(x(1)^2) );
ocp.st( t0 == 0 );
ocp.st( tf == 300 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [   0; 1] );
ocp.st(  x(tf) == [   0; 0] );
ocp.st( -1e-2 <= u <= 1e-2  );
sol = ocp.solve('intervals', 1000, 'degree', 2);

figure(1); 
subplot(311); hold on;
sol.plot(t, x(1))
subplot(312); hold on;
sol.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)