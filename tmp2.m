import yop.*

[t, t0, tf] = independent('t');
s = state('s');
v = state('v');
a = control('a');

% x = [s; v]
% der(x) == [v; a], ...

x = state('x', 2);
s = x(1);
v = x(2);

ocp = yop.ocp();
ocp.min( 0.5 * int(a^2) );
ocp.st(...
    0 <= t <= 1, ... Ska detta tolkas som {t0 >= 0, tf <= 1, tf >= t0} eller {t = [0, 1]} ?
    t0 == 0, tf == 1, ...
    der(s) == v, ...
    der(v) == a, ...    
    s(t0) == s(tf) == 0, ...
    v(t0) == -v(tf) == 1, ...
    s <= 1/9 ...
    );

%% ocp.parse_constraints();
ocp.parse_variables()
ocp.parse_constraints()

%% Extract indices from subsasgn expressions

% Idén är att köra uttrycket med subindices med variables värde satt till
% en uppräkning av dess element, man får sedan ut en vector med de
% relevanta indexen när uttrycket evaluerats. Det kan sedan användas för
% hämta ut rätt subindices.

x = yop.state('x', 10);
x37 = x(3:7);
x46 = x37(2:4);
bc = x26 <= 2;

x.m_value = 1:size(x,1);
indices = evaluate(bc.lhs)

bnd = zeros(size(x));
bnd(indices) = bc.rhs

%%
import yop.*
[t, t0, tf] = independent('t');
x = state('x', 3);
u = control('u');

[bc, y] = rocket_model(x, u);

m0 = 215;
mf = 68;

ocp = yop.ocp();
ocp.max(y.rocket.height(tf));
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