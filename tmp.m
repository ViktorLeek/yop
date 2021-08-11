t = yop.ast_variable();
v1 = yop.ast_variable([5, 5]);
v2 = yop.ast_variable([5, 5]);
v3 = yop.ast_variable();

e1 = (v1-v2)*v3;
e2 = (v1 - v2 / v3).*v1;
rel = e1 == e2;

% print(e1(2,1))
% e1(1) = 1;