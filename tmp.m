t = yop.ast_variable();
v = yop.ast_variable([1, 10]);
v1 = yop.ast_variable([5, 5]);
v2 = yop.ast_variable([5, 5]);
v3 = yop.ast_variable();
v4 = yop.ast_variable([3, 1]);
v5 = yop.ast_variable([3, 1]);

e1 = (v1-v2)*v3;
e2 = (v1 - v2 / v3).*v1;
rel = e1 == e2;

% print(e1(2,1))
% e1(1) = 1;