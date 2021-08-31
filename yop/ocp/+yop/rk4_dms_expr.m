function rk4 = rk4_dms_expr(f, nx, nu, np, nw, steps)
% YOP.RK4_DMS_EXPR - Parameterizes a rk4 integrator that is tailored to 
% integrating mixed timepoints and timedependent values for the 
% direct multiple shooting method.
%   rk4 = rk4integrator(f, T, sx, su, sp)
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
%     rk4 - casadi.Function object. Parameter list rk4(t0, tf, x0, u, p).
%           Integrates the function f from t0, x0 


if nargin == 4
    steps = yop.defaults.rk4_steps;
end

t0  = casadi.MX.sym('t0');
tf  = casadi.MX.sym('tf');
h = (tf-t0)/steps;
x0 = casadi.MX.sym('x0', nx);
u  = casadi.MX.sym('U', nu);
p  = casadi.MX.sym('P', np);
w = casadi.MX.sym('w', nw);
t = t0;
x = x0;
for j=1:steps
    k1 = f(t      , x           , u, p, w);
    k2 = f(t + h/2, x + h/2 * k1, u, p, w);
    k3 = f(t + h/2, x + h/2 * k2, u, p, w);
    k4 = f(t + h  , x + h   * k3, u, p, w);
    t = t + h;
    x = x + h/6*(k1 +2*k2 +2*k3 +k4);
end
rk4 = casadi.Function('F', {t0, tf, x0, u, p, w}, {x});

end