function g = discretize_ode(ocp, N, tau, dt, t, x, u, p)
g = []; % Equality constraints from discretization

% Dynamics
for n=1:N
    dx = x(n).differentiate();
    for r=tau(2:end)
        dxr = dx.evaluate(r);
        tt = t(n).evaluate(r);
        xx = x(n).evaluate(r);
        uu = u(n).evaluate(r);
        pp = p;
        f = ocp.ode.fn(tt, xx, uu, pp);
        g = [g; (dxr - dt*f)];
    end
end

% Continuity
for n=1:N
    gn = x(n).evaluate(1) - x(n+1).evaluate(0);
    g = [g; gn];
end
end