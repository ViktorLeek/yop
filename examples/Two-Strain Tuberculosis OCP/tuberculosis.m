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