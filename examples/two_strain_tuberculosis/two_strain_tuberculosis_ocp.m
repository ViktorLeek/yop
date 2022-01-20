yops Times: t t0 tf
yops State: x size: [6,1] weight: [1e4,1e4,1e4,1e3,1e3,1e3]
yops Ctrls: u size: [2,1]

Npop = 30000;
x0   = [76; 1; 36; 2; 4; 1]*Npop/120;

R = diag([25; 250]);
ocp = yop.ocp('Two-Strain Tuberculosis');
ocp.min( 1e-3*int(x(4) + x(6) + u'*R*u) );
ocp.st( t0==0, tf==5 );
ocp.st( der(x) == tuberculosis(x, u) );
ocp.st(  x(t0) == x0 );
ocp.st( 0.00 <= x <= 30e3 );
ocp.st( 0.05 <= u <= 0.95 );
ocp.st( sum(x) == Npop );

ig = yop.guess(t0, 0, tf, 5, x, x0, u, [0.95; 0.95]);

sol = ocp.solve('intervals', 100, 'guess', ig);

figure(1)
subplot(321); hold on
sol.plot(t, x(1))
subplot(322); hold on
sol.plot(t, x(2))
subplot(323); hold on
sol.plot(t, x(3))
subplot(324); hold on
sol.plot(t, x(4))
subplot(325); hold on
sol.plot(t, x(5))
subplot(326); hold on
sol.plot(t, x(6))

figure(2)
subplot(211); hold on
sol.plot(t, u(1))
subplot(212); hold on
sol.plot(t, u(2))

function dx = tuberculosis(x, u)

S=x(1); T=x(2); L1=x(3); L2=x(4); I1=x(5); I2=x(6); 
u1=u(1); u2=u(2);

beta1 = 13;
beta2 = 13;
mu = 0.0143;
d1 = 0;
d2 = 0;
k1 = 0.5;
k2 = 1;
r1 = 2;
r2 = 1;
p = 0.4;
q = 0.1;
Npop = 30000;
betas = 0.029;
lam = mu*Npop;
m0 = 1;
dm = 0.0749;

dS = lam - (beta1*S*I1 + betas*S*I2)/Npop - mu*S;
dT = u1*r1*L1 - mu*T + (1 - (1-u2)*(p+q))*r2*I1 - (beta2*T*I1 + betas*T*I2)/Npop;
dL1 = (beta1*S*I1 + beta2*T*I1 - betas*L1*I2)/Npop - (mu+k1)*L1 - u1*r1*L1 + (1-u2)*p*r2*I1;
dL2 = (1-u2)*q*r2*I1 - (mu+k2)*L2 + betas*(S+L1+T)*I2/Npop;
dI1 = k1*L1 - (mu+d1)*I1 - r2*I1;
dI2 = k2*L2 - (mu+d2)*I2;
dx  = [dS; dT; dL1; dL2; dI1; dI2];
end