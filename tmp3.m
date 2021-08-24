import yop.*

[t, t0, tf] = independent('t');
x = state('x', 1, 2);
u = control('u');
z = algebraic('z', 1, 4);

eI2_abs = [x(1), x(2), x, u];
eI2_abs = [eI2_abs(4), eI2_abs(3), eI2_abs(1), eI2_abs(2), eI2_abs(4)];
[b, v] = isa_variable([eI2_abs(2:end), z(2:3)+2])
%%
eI2_abs = [true, false, true, true, false];
all([true, true, false])

x = yop.state('x', 3, 3);
eI2_abs = x(1:2, 1:2);

sz = size(eI2_abs.node);
eI2_abs.node.m_value = reshape(1:prod(sz), sz);
forward(eI2_abs)


%%
x = yop.state('x', 3);
y = yop.state('y', 3);

x1 = x(1);
x2 = x(2);
y3 = y(3);

t1 = [x(end-1:end); y];
t2 = [t1(end-2:end); t1(1:end-3)];
bc = t2 <= [1;2;3;4;5];

bd1 = -inf(size(x));
bd2 = -inf(size(x));


% 1) Indexera alla variablers element entydigt

vars = yop.get_variables(bc);

idx0 = 1;
vk = vars{1};
sz = size(vk);
n1 = prod(sz);
I1 = idx0:(idx0+n1-1);
vk.m_value = reshape(I1, sz);
idx0 = idx0 + n1;

vk = vars{2};
sz = size(vk);
n2 = prod(sz);
I2 = idx0:(idx0+n2-1);
vk.m_value = reshape(I2, sz);
idx0 = idx0 + n2;

% 2) Utvärdera uttrycket för att propagera indexen
e = forward_evaluate(bc.lhs);

% 3) För de element som inte är variabler, sätt deras index till något
% värde som uppräkningen inte kan anta
e(~isa_variable(bc.lhs)) = -1;

% Separera indexen på variablerna för att sedan kunna bygga en vektor med
% de relevanta elementen ur högerleden på bivillkoren för den givna
% variabeln
eI1 = e(e >= I1(1) & e <= I1(end));
eI2 = e(e >= I2(1) & e <= I2(end));

% Plocka ut de relevanta indexen för variabeln i den stora vektorn, dvs
% absolut adress i vectorn med alla index
eI1_abs = false(size(e));
for k=1:length(eI1)
    I = find(e==eI1(k));
    eI1_abs(I) = true;
end

eI2_abs = false(size(e));
for k=1:length(eI2)
    I = find(e==eI2(k));
    eI2_abs(I) = true;
end

% Beräkna relative adress i den variabel-lokala vektor
e1_rel = eI1 - I1(1) + 1;
e2_rel = eI2 - I2(1) + 1;

% Sätt variabelns gränser
bd1(e1_rel) = bc.rhs(eI1_abs)
bd2(e2_rel) = bc.rhs(eI2_abs)
