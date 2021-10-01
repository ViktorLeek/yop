function res = load(filename, t, x, u, p)
z = [];

data = load(filename);
sol.x = data.w;

nlp = struct;
nlp.N = data.N;
nlp.d = data.d;
nlp.cp = data.cp;
nlp.tau = [0, casadi.collocation_points(nlp.d, nlp.cp)];
nlp.t0 = yop.cx('t0');
nlp.tf = yop.cx('tf');
nlp.t = yop.collocated_time(nlp.t0, nlp.tf, data.N);
nlp.x = yop.collocated_state(data.nx, nlp.N, nlp.tau);
nlp.z = [];
nlp.u = yop.collocated_control(data.nu, nlp.N);
nlp.p = yop.cx('p', data.np);
nlp.dt = (nlp.tf - nlp.t0)/nlp.N;
nlp.w = vertcat(nlp.t0, nlp.tf, vec(nlp.x), vec(nlp.u), nlp.p);


tt = yop.ocp_var(t);

xx = yop.ocp_var.empty(1,0);
vars = yop.ocp.find_special_nodes(x);
for k=1:length(vars)
    xx(end+1) = yop.ocp_var(vars{k});
end

zz = yop.ocp_var.empty(1,0);
vars = yop.ocp.find_special_nodes(z);
for k=1:length(vars)
    zz(end+1) = yop.ocp_var(vars{k});
end
nlp.tau = [0, casadi.collocation_points(data.d, data.cp)];
uu = yop.ocp_var.empty(1,0);
vars = yop.ocp.find_special_nodes(u);
for k=1:length(vars)
    uu(end+1) = yop.ocp_var(vars{k});
end

pp = yop.ocp_var.empty(1,0);
vars = yop.ocp.find_special_nodes(p);
for k=1:length(vars)
    pp(end+1) = yop.ocp_var(vars{k});
end

res = yop.ocp_sol(sol, [], nlp, tt, xx, zz, uu, pp);
end