yopvar times: t t0 tf states: p q controls: u

zeta = 0.084; % per day
b = 5.85;     % per day
d = 0.00873;  % per mm^2 per day
G = 0.15;     % per mg of dose per day
mu = 0.02;    % per day
a = 75; 
A = 15;

x_max = ((b-mu)/d)^(3/2);
x_min = 0.1;
u_max = a;
u_min = 0;

ocp = yop.ocp('Tumor Antiangiogenesis');
ocp.min(p(tf));
ocp.st( ...
    0.1 <= tf <= 5, ...
    der(p) == -zeta*p*log(p/q), ...
    der(q) == q * (b - (mu + (d*(p^(2/3))) + G*u)), ...
    p(t0) == 0.50*x_max, ...
    q(t0) == 0.25*x_max, ...
    x_min <= [p; q] <= x_max, ...
    u_min <=    u   <= u_max, ...
    int(u) <= A ...
    );

sol = ocp.solve('intervals', 150);

figure(1)
subplot(311); hold on;
sol.plot(t, p);
subplot(312); hold on;
sol.plot(t, q);
subplot(313); hold on;
sol.stairs(t, u);