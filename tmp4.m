syms v1 v2 v3 v4 x1 x2
vv = [v1; v2; v3; v4];

v = yop.state('v', 4);
v.m_value = vv;

x = yop.state('x', 2);
x.m_value = [x1; x2];

[t, t0, tf] = yop.independent();


v12 = v(1:2);
r1 = 2*[x(1); x.^2; v12(1)*2; v12(2)] <= [x(1); x(t0); v12(1)*2; v12(2)]; 

tmp = split_vars_and_exprs(r1);

clc
disp('Constraint: ')
disp(forward_evaluate(r1));
disp('is split into: ')


for k=1:length(tmp)
    disp(['[k=' num2str(k) ']'])
    disp(forward_evaluate(tmp{k}))
end


