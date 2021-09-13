%% Formulation 3
% Time
[t, t0, tf] = yop.time('t');

% Rocket parameters
D0 = 0.01227; beta = 0.145e-3; c = 2060;

% Constants
g0 = 9.81; r0 = 6.371e6;

% Rocket model
h  = yop.state('h');  % Rocket height
v  = yop.state('v');  % Rocket speed
m  = yop.state('m');  % Rocket mass
Wf = yop.control('Wf');  % Rocket fuel massflow

% Drag force and gravitational acceleration
F_D = D0 * exp(-beta*h) * v^2;
g   = g0*(r0/(r0+h))^2;

% Mass boundaries
m_min = 68; m_max = 215; 

% Control boundaries
Wfmin = 0; Wfmax = 9.5;

% Optimal control problem
grp = yop.ocp();
grp.max( h(tf) );
grp.st( ...
    ... Rocket model
    der(h) == v              , ...
    der(v) == (Wf*c-F_D)/m-g , ...
    der(m) == -Wf            , ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_max, ...
    ... Box constraints
    h >= 0, ...
    v >=0, ...
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

%%
t.m_value = sym('t');
t0.m_value = sym('t0');
tf.m_value = sym('tf');
h.m_value = sym('h');
v.m_value = sym('v');
m.m_value = sym('m');
Wf.m_value = sym('W_f');

e1 = (Wf*c-F_D)/m-g;
e1(1) = v;

% cc = { ...
%     [m_min; Wfmin] <= [m; Wf] <= [m_max; Wfmax], ...
%     [v; der(v)] >= [0; 10], ...
%     e1 <= 100, ...
%     der(v) == (Wf*c-F_D)/m-g ...
%     };

cc = { ...
    [m_min; Wfmin] <= [m(t0); Wf(tf)] <= [m_max; Wfmax], ...
    [v(t==2); der(v)] >= [0; 10], ...
    e1 <= 100, der(v) == (Wf*c-F_D)/m-g ...
    };
srf = yop.ocp.to_srf(cc);
[bc, gc] = yop.ocp.get_box(srf);
cbox = yop.ocp.canonicalize(bc)
ubox = yop.ocp.unique_box(cbox)
[bc_t, bc_t0, bc_tf] = yop.ocp.timed_box(ubox)

%%
tmp = cbox;
for k=1:length(tmp)
    %     draw(tmp{k}.lhs)
    %     disp('----------')
    propagate_value(tmp{k})
end
