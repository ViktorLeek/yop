%% LD2
load('models/liu_diesel_2/params.mat')

% simulate
u_num = [15; 1; 0];
f = @(t, x) liu_diesel_2(x, u_num, 1200, ice_param);
[t_sim, x_sim] = ode15s(f, [0, 10], [2e5; 2e5; 2e5; 9000]);
x_num = x_sim(end,:)';

% evaluate derivative at end of interval
x = yop.state('x', 4, 1);
u = yop.control('u', 3, 1);
disp('---- exact ----')
dx_num = liu_diesel_2(x_num, u_num, 1200, ice_param)

% Test evaluation mode to see that values match up
[dx, y, h] = liu_diesel_2(x, u, 1200, ice_param);
x.m_value = x_num;
u.m_value = u_num;

disp('---- forward ----')
tic()
[topsort, ~, n_elem] = topological_sort(dx);

for k=1:(n_elem-1)
    forward(topsort{k});
end
forward(topsort{n_elem}) % returns value
toc()

disp('---- recursive ----')
tic()
evaluate(dx)
toc()

%% Nested calls

x0 = x;
dx0 = liu_diesel_2(x0, u, 1200, ice_param);
x1 = x0 + 1e-5*dx0;
dx1 = liu_diesel_2(x1, u, 1200, ice_param);
x2 = x1 + 1e-5*dx1;
dx2 = liu_diesel_2(x2, u, 1200, ice_param);
x3 = x2 + 1e-5*dx2;
dx3 = liu_diesel_2(x3, u, 1200, ice_param);
x4 = x3 + 1e-5*dx3;
dx4 = liu_diesel_2(x4, u, 1200, ice_param);
x5 = x4 + 1e-5*dx4;
dx5 = liu_diesel_2(x5, u, 1200, ice_param);
x6 = x5 + 1e-5*dx5;
dx6 = liu_diesel_2(x6, u, 1200, ice_param);

clc
disp('---- forward ----')

tic()
[topsort, ~, n_elem] = topological_sort(dx6);
% toc()
% n_elem

% tic()
for k=1:(n_elem-1)
    forward(topsort{k});
end
forward(topsort{n_elem});
toc()

% disp('---- recursive ----')
% tic()
% evaluate(dx2)
% toc()

%% LDE
x = yop.ast_variable([5, 1]);
u = yop.ast_variable([3, 1]);

x0 = [83.7758; 1.0143e+05; 1.0975e+05; 2.0502e+03; 0];
u_num = [15; 0; 0];
gensetModel(x0, u_num)

dx = gensetModel(x, u);
x.value = x0;
u.value = u_num;
evaluate(dx)

%% Goddard's rocket problem
x = yop.ast_variable([3, 1]);
u = yop.ast_variable();

x0 = [100; 1000; 200];
u_num = 8;
rocket_model(x0, u_num)
dx = rocket_model(x, u);
x.value = x0;
u.value = u_num;
evaluate(dx)

%%
syms xx1 xx2
x1 = yop.state('x1');
x2 = yop.state('x2');
expr = abs(x1) + x2;
[topsort, visited] = topological_sort(expr)
x1.value = xx1;
x2.value = xx2;