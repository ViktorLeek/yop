function rk4 = rk4_integrator(f, nx, nu, np, steps)
% YOP.RK4_INTEGRATOR - Parameterizes a rk4 integrator.
%   
%   Parameters:
%         f - ode rhs, a function handle or casadi.Function object, with
%             parameterlist f(t,x,u,p)
%        sx - Size of state vector [rows, cols]
%        su - Size of control vector [rows, cols]
%        sp - Size of parameter vector [rows, cols]
% 
%   (Optional)
%     steps - The number of Runge-Kutta steps to take, i.e. the number of
%             times to apply RK4. The more steps, the more the local error
%             is reduced, bit it is also more expensive.
% 
%   Return value:
%     rk4 - casadi.Function object. Parameter list rk4(t0, tf, p0, pf, x0, u, p).
%           Integrates the function f from t0, x0 


if nargin == 4
    steps = yop.defaults.rk4_steps;
end

p0  = casadi.MX.sym('t0'); % Parameter t0 for entire problem
pf  = casadi.MX.sym('tf'); % Parameter tf for entire problem
t0  = casadi.MX.sym('t0'); % Start of integration
tf  = casadi.MX.sym('tf'); % End of integration
x0 = casadi.MX.sym('x0', nx);
u  = casadi.MX.sym('U', nu);
p  = casadi.MX.sym('P', np);
t = t0;
x = x0;
h = (tf-t0)/steps;
for j=1:steps
    k1 = f(p0, pf, t      , x           , u, p);
    k2 = f(p0, pf, t + h/2, x + h/2 * k1, u, p);
    k3 = f(p0, pf, t + h/2, x + h/2 * k2, u, p);
    k4 = f(p0, pf, t + h  , x + h   * k3, u, p);
    t = t + h;
    x = x + h/6*(k1 +2*k2 +2*k3 +k4);
end
rk4 = casadi.Function('F', {t0, tf, p0, pf, x0, u, p}, {x});

end