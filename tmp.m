%% LD2
load('models/ld2/params.mat')
x = yop.ast_variable([4, 1]);
u = yop.ast_variable([3, 1]);

x0 = [2e5; 2e5; 2e5; 9000]; 
u0 = [15; 1; 0];
f = @(t, x) liu_diesel_2(x, u0, 1200, ice_param);
[t_sim, x_sim] = ode15s(f, [0, 10], x0);
x_num = x_sim(end,:)';
dx_num = liu_diesel_2(x_num, u0, 1200, ice_param)

[dX, y, h] = liu_diesel_2(x, u, 1200, ice_param);
x.value = x_num;
u.value = u0;
evaluate(dX)

%% LDE
x = yop.ast_variable([5, 1]);
u = yop.ast_variable([3, 1]);

x0 = [83.7758; 1.0143e+05; 1.0975e+05; 2.0502e+03; 0];
u0 = [15; 0; 0];
gensetModel(x0, u0)

dx = gensetModel(x, u);
x.value = x0;
u.value = u0;
evaluate(dx)

%% Goddard's rocket problem
x = yop.ast_variable([3, 1]);
u = yop.ast_variable();

x0 = [100; 1000; 200];
u0 = 8;
rocket_model(x0, u0)
dx = rocket_model(x, u);
x.value = x0;
u.value = u0;
evaluate(dx)