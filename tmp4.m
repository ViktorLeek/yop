syms v1 v2 v3 v4 x1 x2 x3 x4
vv = [v1; v2; v3; v4];

v = yop.state('v', 4);
v.m_value = vv;

x = yop.state('x', 4);
x.m_value = [x1; x2; x3; x4];

[t, t0, tf] = yop.independent();

mix = [v(1); x(1); v(2); x(2); v(3); x(3); v(4); x(4)];

constraint = 1 <= [mix; mix(3)*2] <= 1; 

srf = to_srf(constraint);
srf_vne = {};
for k=1:length(srf)
    vne = split_vars_and_exprs(srf{k});
    srf_vne = {srf_vne{:}, vne{:}};
end

tmp = {};
for k=1:length(srf_vne)
    tk = split_vars(srf_vne{k});
    tmp = {tmp{:}, tk{:}};
end

clc
disp('Constraint: ')
disp(forward_evaluate(constraint));
disp('is transformed into: ')


for k=1:length(tmp)
    disp(['[k=' num2str(k) ']'])
    disp(forward_evaluate(tmp{k}))
end
