% Debugging using syms
syms tt xx1 xx2 xx3 uu
xx = [xx1; xx2; xx3];

import yop.*
t = ast_variable();
x = ast_variable([3, 1]);
u = ast_variable();
[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp();
ocp.max(rocket.height);
ocp.st(...
    der(x) == rocket_model(x, u)      , ...
    ... t(0)  == 0                     , ...
    rocket.height   >= 0              , ...
    rocket.velocity >= 0              , ...
    68 <=     rocket.mass     <= 215  , ...
     0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

x.value = xx;
u.value = uu;