load('models/ld2/params.mat')
t = yop.ast_variable();
x = yop.ast_variable([4, 1]);
u = yop.ast_variable([3, 1]);

[dX, signals, constraints, param] = liu_diesel_2(x, u, 1200, ice_param);
