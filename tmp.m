v1 = yop.variable([5, 1]);
v2 = yop.variable([5, 1]);
v3 = yop.variable();

e1 = (v1-v2)*v3;
e2 = (v1 - v2 / v3).*v1;
rel = e1 == e2;
print(e1)