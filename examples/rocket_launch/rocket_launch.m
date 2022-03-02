%% Variables and parameters
yops Times: t t0 tf 
yops States: rx ry rz vx vy vz mf nominal: [1e6,1e6,1e6,1e3,1e3,1e3,1e5]
yops Controls: u size: [3,1] int: 1
r=[rx;ry;rz]; v=[vx;vy;vz]; x=[r;v;mf];
[srb,e1,e2,Re,w,mu] = problem_data();

%% Initial value
r0  = [Re*cosd(28.5); 0; Re*sind(28.5)];
v0  = -cross(w,r0) + 5*r0/norm(r0);
mf0 = 9*srb.mp + e1.mp + e2.mp;
x0  = [r0; v0; mf0];

%% Final value
d2r = pi/180;
final_orbit = [24361140; 0.7308; 28.5*d2r; 269.8*d2r; 130.5*d2r; 0];

%% Initial guess
ivp1 = yop.ivp();
ivp1.add(t0==0, tf==srb.tb);
ivp1.add( der(x) == vehicle1(x,u) );
ivp1.add( x(t0) == x0 );
ivp1.add( u == [1;0;0] );
sim1 = ivp1.solve();

ivp2 = yop.ivp();
ivp2.add(t0==srb.tb, tf==2*srb.tb);
ivp2.add( der(x) == vehicle2(x,u) );
ivp2.add( x(t0) == sim1.value(x(tf)) );
ivp2.add( u == [1;0;0] );
sim2 = ivp2.solve();

ivp3 = yop.ivp();
ivp3.add(t0==2*srb.tb, tf==e1.tb);
ivp3.add( der(x) == vehicle3(x,u) );
ivp3.add( x(t0) == sim2.value(x(tf)) );
ivp3.add( u == [0;1;0] );
sim3 = ivp3.solve();

ivp4 = yop.ivp();
ivp4.add(t0==e1.tb, tf==e1.tb+e2.tb);
ivp4.add( der(x) == vehicle4(x,u) );
ivp4.add( x(t0) == sim3.value(x(tf)) );
ivp4.add( u == [0;1;0] );
sim4 = ivp4.solve();

%% OCP

p1 = yop.ocp();
p1.st( 0 == t0 < tf == srb.tb );
p1.st( der(x) == vehicle1(x,u) );
p1.st( x(t0) == x0 );
p1.st( norm(u)^2 == 1 );
p1.st( norm(r)^2 >= Re^2 );
p1.guess = sim1;

p2 = yop.ocp();
p2.st( srb.tb == t0 < tf == 2*srb.tb );
p2.st( der(x) == vehicle2(x,u) );
p2.st( norm(u)^2 == 1 );
p2.st( norm(r)^2 >= Re^2 );
p2.guess = sim2;

p3 = yop.ocp();
p3.st( 2*srb.tb == t0 < tf == e1.tb );
p3.st( der(x) == vehicle3(x,u) );
p3.st( norm(u)^2 == 1 );
p3.st( norm(r)^2 >= Re^2 );
p3.guess = sim3;

p4 = yop.ocp();
p4.max( mf(tf) );
p4.max( norm(r(tf))^2 );
p4.st( e1.tb == t0 < tf <= e1.tb+e2.tb );
p4.st( der(x) == vehicle4(x,u) );
p4.st( norm(u)^2 == 1 );
p4.st( norm(r)^2 >= Re^2 );
p4.st( state2orbit(x(tf),mu) == final_orbit );
p4.guess = sim4;

%%
% ocp = p1 + p2 + p3 + p4;
% sol = ocp.solve('ival',50);
% vehicle1(x0,[0.4;0.6;0])

function dx = vehicle1(x, u)
% Launch vehicle dynamics for phase 1

r = x(1:3);
v = x(4:6);
mf = x(7);

[srb,e1,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data();

vr  = v - cross(w,r); % Earth relative velocity
% vr  = v + cross(w,r); % Earth relative velocity
h   = max(norm(r) - Re, 0);
rho = rho0*exp(-h/h0);
D   = 0.5 * Af * cd * rho * norm(vr) * vr;
T   = 6*srb.T + e1.T;
m   = 9*srb.m0 + e1.m0 + e2.m0 + mf + mpl;
dm  = -6*srb.T/(g0*srb.I) - e1.T/(g0*e1.I);
dv  = T/m*u - D/m - mu/norm(r)^3 * r;
dx  = [v; dv; dm];

end

function dx = vehicle2(x, u)
% Launch vehicle dynamics for phase 2

r = x(1:3);
v = x(4:6);
mf = x(7);

[srb,e1,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data();

vr  = v - cross(w,r); % Earth relative velocity
% vr  = v + cross(w,r); % Earth relative velocity
h   = max(norm(r) - Re, 0);
rho = rho0*exp(-h/h0);
D   = 0.5 * Af * cd * rho * norm(vr) * vr;
T   = 3*srb.T + e1.T;
m   = 3*srb.m0 + e1.m0 + e2.m0 + mf + mpl;
dm  = -3*srb.T/(g0*srb.I) - e1.T/(g0*e1.I);
dv  = T/m*u - D/m - mu/norm(r)^3 * r;
dx  = [v; dv; dm];

end

function dx = vehicle3(x, u)
r = x(1:3);
v = x(4:6);
mf = x(7);

[~,e1,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data();

vr  = v - cross(w,r); % Earth relative velocity
% vr  = v + cross(w,r); % Earth relative velocity
h   = max(norm(r) - Re, 0);
rho = rho0*exp(-h/h0);
D   = 0.5 * Af * cd * rho * norm(vr) * vr;
T   = e1.T;
m   = e1.m0 + e2.m0 + mf + mpl;
dm  = -e1.T/(g0*e1.I);
dv  = T/m*u - D/m - mu/norm(r)^3 * r;
dx  = [v; dv; dm];

end

function dx = vehicle4(x, u)
r = x(1:3);
v = x(4:6);
mf = x(7);

[~,~,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data();

vr  = v - cross(w,r); % Earth relative velocity
% vr  = v + cross(w,r); % Earth relative velocity
h   = max(norm(r) - Re, 0);
rho = rho0*exp(-h/h0);
D   = 0.5 * Af * cd * rho * norm(vr) * vr;
T  = e2.T;
m  = e2.m0 + mf + mpl;
dm = -e2.T/(g0*e2.I);
dv = T/m*u - D/m - mu/norm(r)^3 * r;
dx = [v; dv; dm];

end


function dx = vehicle(x, u, phase)
r = x(1:3);
v = x(4:6);
mf = x(7);

[srb,e1,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data();

% omegaMatrix = w(3)*[0 -1 0;1 0 0;0 0 0];
% omegacrossr = transpose(r'*omegaMatrix.');
% vr = v-omegacrossr;

vr  = v - cross(w,r); % Earth relative velocity
h   = norm(r) - Re;
rho = rho0*exp(-h/h0);
D   = 0.5 * Af * cd * rho * norm(vr) * vr;
switch phase
    case 1
        T  = 6*srb.T + e1.T;
        m  = 9*srb.m0 + e1.m0 + e2.m0 + mf + mpl;
        dm = -6*srb.T/(g0*srb.I) - e1.T/(g0*e1.I);
    case 2
        T  = 3*srb.T + e1.T;
        m  = 3*srb.m0 + e1.m0 + e2.m0 + mf + mpl;
        dm = -3*srb.T/(g0*srb.I) - e1.T/(g0*e1.I);
    case 3
        T  = e1.T;
        m  = e1.m0 + e2.m0 + mf + mpl;
        dm = -e1.T/(g0*e1.I);
    case 4
        T  = e2.T;
        m  = e2.m0 + mf + mpl;
        dm = -e2.T/(g0*e2.I);
end
dv = T/m*u - D/m - mu/norm(r)^3 * r;
dx = [v; dv; dm];

end

function [srb,e1,e2,Re,w,mu,g0,rho0,h0,Af,cd,mpl] = problem_data()

g0 = 9.80665;
w = [0; 0; 7.29211585e-5]; % Angular velocity of earth
Re = 6378145; % Earth radius
rho0 = 1.225; % Sea level density
h0 = 7200; 
Af = 4*pi; % Frontal area
cd = 0.5; % Drag coefficient
mu = 3.986012e14;
mpl = 4164; % payload

srb.m  = 19290;  % SRB Total mass
srb.mp = 17010;  % SRB propellant mass
srb.m0 = srb.m - srb.mp; % SRB dry mass
srb.T  = 628500; % SRB Thrust 
srb.tb = 75.2;   % SRB Burn time
srb.I  = srb.T/(g0 * srb.mp/srb.tb); % SRB Specific impulse

e1.m  = 104380;  % Main engine total mass
e1.mp = 95550;   % Main engine propellant mass
e1.m0 = e1.m - e1.mp; % Main engine dry mass
e1.T  = 1083100; % Main engine thrust
e1.tb = 261;     % Main engine burn time
e1.I  = e1.T/(g0 * e1.mp/e1.tb); % Main engine specific impulse

e2.m  = 19300;  % Secondary engine total mass
e2.mp = 16820;  % Secondary engine propellant mass
e2.m0 = e2.m - e2.mp; % Secondary engine dry weight
e2.T  = 110094; % Secondary engine thrust
e2.tb = 700;    % Secondary engine burn time
e2.I  = e2.T/(g0 * e2.mp/e2.tb); % Secondary engine specific impulse

end

function rv = orbit2rv(oe,mu)
a=oe(1); e=oe(2); i=oe(3); O=oe(4); w=oe(5); n=oe(6); s=@sin; c=@cos;
pp  = a*(1-e^2);
rr  = pp/(1+e*c(n));
rv = [rr*c(n); rr*s(n); 0];
vv = sqrt(mu/pp)*[-s(n); e+c(n); 0];
RR = [c(O)*c(w)-s(O)*s(w)*c(i), -c(O)*s(w)-s(O)*c(w)*c(i),  s(O)*s(i);
      s(O)*c(w)+c(O)*s(w)*c(i), -s(O)*s(w)+c(O)*c(w)*c(i), -c(O)*s(i);
      s(w)*s(i)               ,  c(w)*s(i)               ,  c(i)    ];
r = RR*rv;
v = RR*vv;
rv = [r;v];
end

function oe = state2orbit(x,mu)
rv = x(1:3);
vv = x(4:6);
K  = [0;0;1];
hv = cross(rv, vv);
nv = cross(K, hv);
n  = sqrt(nv'*nv);
h2 = hv'*hv;
v2 = vv'*vv;
r  = sqrt(rv'*rv);
ev = 1/mu *( (v2-mu/r)*rv - (rv'*vv)*vv ); 
p  = h2/mu;
e  = sqrt(ev'*ev);
a  = if_else(e==1, p/0.001, p/(1-e*e));
i  = if_else(h2==0, acos(hv(3)/0.001) ,acos(hv(3)/sqrt(h2)));
Om = if_else(n==0, acos(nv(1)/0.001), acos(nv(1)/n));
om = if_else(n==0|e==0, acos(nv.'*ev/0.001), acos(nv.'*ev/n/e)); 
nu = if_else(e==0|r==0, acos(ev.'*rv/0.001), acos(ev.'*rv/e/r)); 
oe = [a; e; i; Om; om; nu];
end