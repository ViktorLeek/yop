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

disp('---- forward ----')
tic()
[topsort, ~, n_elem] = topological_sort(dx4);
n_elem

for k=1:n_elem
    forward(topsort{k});
end

topsort{n_elem}.m_value
toc()