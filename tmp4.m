syms v1 v2 v3 v4 x1 x2 x3 x4
vv = [v1; v2; v3; v4];

v = yop.state('v', 4);
v.m_value = vv;

x = yop.state('x', 4);
x.m_value = [x1; x2; x3; x4];

[t, t0, tf] = yop.independent();

v33 = v(3);
v11 = v(1);
x44 = x(4);

mix = [v(1); x(1); v(2); x(2); v33(t0); x(3); v(4); x(4)];
mix(1) = v11(t==4);
mix(end) = x44(tf);
m67 = mix([6,7]);
c1 = 1 <= [mix; mix(3)*2; exp(norm(m67)); mix] <= [v(2); v(3).^x(2); 2*ones(16,1)]; 
% 
% c1 = [x(1); x(2); v3(t0); v(2)*x(1)] <= 1;

constraints = {c1};
srf = yop.to_srf(constraints);
hsrf = yop.to_hsrf(srf.get_relations);
vnf = yop.to_vnf(hsrf);
dtp = yop.to_dtp(vnf);

clc
disp('Constraint:')
disp(forward_evaluate(c1));
disp('is transformed into:')

disp('Box constraints:')

for k=1:length(dtp.vn_t)
    disp(evaluate(dtp.vn_t{k}))
end

for k=1:length(dtp.nv_t)
    disp(evaluate(dtp.nv_t{k}))
end

disp('Initial box constraints:')

for k=1:length(dtp.vn_t0)
    disp(evaluate(dtp.vn_t0{k}))
end

for k=1:length(dtp.nv_t0)
    disp(evaluate(dtp.nv_t0{k}))
end

disp('Final box constraints:')

for k=1:length(dtp.vn_tf)
    disp(evaluate(dtp.vn_tf{k}))
end

for k=1:length(dtp.nv_tf)
    disp(evaluate(dtp.nv_tf{k}))
end

disp('General constraints:')
for k=1:length(dtp.ve)
    disp(evaluate(dtp.ve{k}))
end

for k=1:length(dtp.ev)
    disp(evaluate(dtp.ev{k}))
end

for k=1:length(dtp.vv)
    disp(evaluate(hsrf.vv{k}))
end

for k=1:length(dtp.ee)
    disp(evaluate(dtp.ee{k}))
end







