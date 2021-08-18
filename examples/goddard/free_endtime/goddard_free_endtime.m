%%
syms tt xx1 xx2 xx3 uu
xx = [xx1; xx2; xx3];

import yop.*
[t, t0, tf] = independent('t');
x = state('x', 3);
u = control('u');

[dx, y] = rocket_model(x, u);

m0 = 215;
mf = 68;

ocp = yop.ocp();
ocp.max(y.rocket.height);
ocp.st(...
    ... dynamics 
    der(x) == rocket_model(x, u), ...
    ... initial condition
    t0  == 0, ...
    y.rocket.height(t0) == 0, ...
    y.rocket.velocity(t0) == 0, ...
    y.rocket.mass(t0) == m0, ...
    ... Box constraints
    y.rocket.height >= 0, ...
    y.rocket.velocity >= 0, ...
    mf <= y.rocket.mass <= m0, ...
     0 <= y.rocket.fuel_mass_flow <= 9.5 ...
    );

x.m_value = xx;
u.m_value = uu;

% rocket.height(t0);
% rocket.height(tf);
% rocket.height(t==t0);
% rocket.height(t==tf);
% rocket.height(t==4);