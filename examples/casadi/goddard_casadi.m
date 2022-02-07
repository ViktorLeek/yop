% Goddard Rocket CasADi implementation.
% Modified from the following example on their homepage.
% https://github.com/casadi/casadi/blob/master/docs/examples/matlab/direct_collocation.m
% https://web.casadi.org/docs/#a-simple-test-problem

import casadi.*

% Degree of interpolating polynomial
d = 5;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

% Time horizon
T = MX.sym('T');

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x = [x1; x2; x3];
u = SX.sym('u');

% Model equations
xdot = rocket_model(x,u);

% Continuous time dynamics
f = Function('f', {x, u}, {xdot});

% Control discretization
N = 50; % number of control intervals
h = T/N;

% Start with an empty NLP
w={}; lbw = []; ubw = []; g={};

% Rocket parameters
m0 = 215;
mf = 68;
Fm = 9.5;

% Final time condition
w = {w{:}, T};
lbw = [lbw; 0];
ubw = [ubw; inf];

% Initial conditions
Xk = MX.sym('X0', 3);
w = {w{:}, Xk};
lbw = [lbw; 0; 0; m0];
ubw = [ubw; 0; 0; m0];

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw; 0];
    ubw = [ubw;  Fm];

    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 3);
        w = {w{:}, Xkj{j}};
        lbw = [lbw;   0;   0; mf];
        ubw = [ubw; inf; inf; m0];
    end

    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xk;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkj{r};
       end

       % Append collocation equations
       g = {g{:}, h*f(Xkj{j},Uk) - xp};

       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};
    end

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 3);
    w = {w{:}, Xk};
    lbw = [lbw;   0;    0; mf];
    ubw = [ubw; inf;  inf; m0];

    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
end

% Create an NLP solver
prob = struct('f', -Xk(2), 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', ones(size(vertcat(w{:}))), 'lbx', lbw, 'ubx', ubw,...
            'lbg', zeros(size(vertcat(g{:}))), ...
            'ubg', zeros(size(vertcat(g{:}))));
w_opt = full(sol.x);

% Plot the solution
T_opt  = w_opt(1);
x1_opt = w_opt(2:4+3*d:end);
x2_opt = w_opt(3:4+3*d:end);
x3_opt = w_opt(4:4+3*d:end);
u_opt  = w_opt(5:4+3*d:end);
tgrid  = linspace(0, T_opt, N+1);

figure(1)
subplot(311); hold on
plot(tgrid, x1_opt)
xlabel('Time'); ylabel('Velocity')

subplot(312); hold on
plot(tgrid, x2_opt)
xlabel('Time'); ylabel('Height')

subplot(313); hold on
plot(tgrid, x3_opt)
xlabel('Time'); ylabel('Mass')

figure(2); hold on
stairs(tgrid, [u_opt; nan])
xlabel('Time'); ylabel('F (Control)')
