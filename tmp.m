load('models/ld2/params.mat')
t = yop.ast_variable();
x = yop.ast_variable([1, 4]);
u = yop.ast_variable([1, 3]);

tic();
[dX, y, h] = liu_diesel_2(x, u, 1200, ice_param);
dX2 = liu_diesel_2(x+dX, u+1, 1200, ice_param);
dX3 = liu_diesel_2(x+dX2, u+1, 1200, ice_param);
dX4 = liu_diesel_2(x+dX3, u+1, 1200, ice_param);
dX5 = liu_diesel_2(x+dX4, u+1, 1200, ice_param);
dX6 = liu_diesel_2(x+dX5, u+1, 1200, ice_param);
dX7 = liu_diesel_2(x+dX6, u+1, 1200, ice_param);
toc();

%%
x0 = [2e5; 2e5; 2e5; 9000]; 
u0 = [15; 1; 0];
f = @(t, x) liu_diesel_2(x, u0, 1200, ice_param);
[t_sim, x_sim] = ode15s(f, [0, 10], x0);
x_num = x_sim(end,:)';
tic
dx_num = liu_diesel_2(x_num, u0, 1200, ice_param);
toc

%%
t = yop.ast_variable();
x = yop.ast_variable([4, 1]);
u = yop.ast_variable([3, 1]);
[dX, y, h] = liu_diesel_2(x, u, 1200, ice_param);

x.value = x_num;
u.value = u0;
tic
evaluate(dX);
toc

%%
[dX, y, h] = liu_diesel_2(x, u, 1200, ice_param);
dX2 = liu_diesel_2(x + 1e-4*dX, u, 1200, ice_param);
dX3 = liu_diesel_2(x + 1e-4*dX2, u, 1200, ice_param);
dX4 = liu_diesel_2(x + 1e-4*dX3, u, 1200, ice_param);
dX5 = liu_diesel_2(x + 1e-4*dX4, u, 1200, ice_param);
dX6 = liu_diesel_2(x + 1e-4*dX5, u, 1200, ice_param);
dX7 = liu_diesel_2(x + 1e-4*dX6, u, 1200, ice_param);
evaluate(dX3)

%%
t = yop.ast_variable();
t.value = 2;
e = t + t;
evaluate(e)






