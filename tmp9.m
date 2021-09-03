[t, t0, tf] = yop.time('t');
p = yop.state('p');   % position
s = yop.state('s');   % speed
a = yop.control('a'); % acceleration
l = 1/9;

x = yop.state('x', 2);
p = x(1);
s = x(2);

d1 = yop.state('d1'); 
d2 = yop.state('d2', 2); 

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( d2(1) + 0.5*int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(s) == a, ... 
    d1(t0) == 0, ...
    der(d2(2)) == 1, ...
    der(p) == s, ...
    p(t0) == 1, d2(tf) == 0, s(t0) ==  1, ...
    p(tf) == 1, s(tf) == -1, ...
    p <= l ...
    );
ocp.build().present();

%% Permute state vector and set rhs
ocp.set_sym;
obj = ocp;

ode_var = obj.odes(1).var;
ode_expr = obj.odes(1).expr;
for k=2:length(obj.odes)
    ode_var = [ode_var(:); obj.odes(k).var(:)];
    ode_expr = [ode_expr(:); obj.odes(k).expr(:)];
end

x_cell = arrayfun(@(e) e.var, obj.states, 'UniformOutput', false);
[re, nr] = yop.reaching_elems(ode_var, x_cell);

% 3 cases: variables that reaches fully, partly, and not at all.
%   1,2) fully and partly
for k=1:length(re)
    if length(re(k).enum) ~= length(re(k).reaching)
        [~, idx] = setdiff(re(k).enum, re(k).reaching);
        vk = re(k).var(idx);
        ode_var = [ode_var(:); vk(:)];
        ode_expr = [ode_expr(:); zeros(size(vk(:)))];
    end
end

%   3) Do not reach at all
for k=1:length(nr)
    vk = nr(k).var;
    ode_var = [ode_var(:); vk(:)];
    ode_expr = [ode_expr(:); zeros(size(vk(:)))];
end

evaluate(ode_var == ode_expr) 

obj.set_mx();
[t0, tf, t, x, u, p] = get_paramlist(obj);
x_fn = casadi.Function('x', {t0, tf, t, x, u, p}, {evaluate(ode_var)});
Px_expr = x_fn(t0, tf, t, x, u, p);
P = casadi.Function('P', {x}, {Px_expr});


%%
states = cell(size(obj.states));
for k=1:length(obj.states)
    states{k} = obj.states(k).var;
end

[~, ode_id] = isa_variable(ode_var);
state_id = yop.get_ids(states);

[diff_id, diff_idx] = setdiff(state_id, ode_id);
x_no_ode = states(diff_idx);
tmp = vertcat(x_no_ode{:});

ode_var = [ode_var(:); tmp(:)];
ode_expr = [ode_expr(:); zeros(size(tmp(:)))];

obj.set_mx();
[t0, tf, t, x, u, p] = get_paramlist(obj);
x_fn = casadi.Function('x', {t0, tf, t, x, u, p}, {evaluate(ode_var)});

Px_expr = x_fn(t0, tf, t, x, u, p);

P = casadi.Function('P', {x}, {Px_expr});

% der() == 0;


%% Process timepoints
%% Process integrals


