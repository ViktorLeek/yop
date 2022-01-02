function nlp = direct_collocation(ocp, N, d, cp)

% If horizon is not fixed, T0 and Tf are nonsense. But since they are only
% to be used for a fixed horizon, that is OK.
% [~, T0, Tf] = ocp.fixed_horizon();

T0=[]; Tf=[]; t=[]; t0=[]; tf=[]; p=[]; dt=[]; tau=[]; tp=[]; I=[]; D=[]; 
g=[]; g_ub=[]; g_lb=[]; w_ub=[]; w_lb=[];
t = yop.interpolating_poly.empty(N+1, 0);
x = yop.interpolating_poly.empty(N+1, 0);
z = yop.interpolating_poly.empty(N  , 0);
u = yop.interpolating_poly.empty(N  , 0);

init_collocation();
init_variables();
special_nodes();
discretize_dynamics();
pointcons();
pathcons();
hard_pathcons();
for pk = ocp.path_ival
    ival_pathcons(pk);
end

w = vertcat(t0, tf, x.vec(), z.vec(), u.vec(), p);
box_bnd();

J = ocp.objective.fn(t0, tf, p, tp, I);

t0fn = casadi.Function('t0', {w}, {t0});
tffn = casadi.Function('tf', {w}, {tf});
tfn  = casadi.Function('t' , {w}, {t.mat()});
xfn  = casadi.Function('x' , {w}, {x.mat()});
zfn  = casadi.Function('z' , {w}, {z.mat()});
ufn  = casadi.Function('u' , {w}, {u.mat()});
pfn  = casadi.Function('p' , {w}, {p});

nlp = struct;
nlp.J = J;
nlp.w = w;
nlp.w_ub = w_ub;
nlp.w_lb = w_lb;
nlp.g = g;
nlp.g_ub = g_ub;
nlp.g_lb = g_lb;
nlp.t0 = @(w) full(t0fn(w));
nlp.tf = @(w) full(tffn(w));
nlp.t  = @(w) full(tfn(w));
nlp.x  = @(w) full(xfn(w));
nlp.z  = @(w) full(zfn(w));
nlp.u  = @(w) full(ufn(w));
nlp.p  = @(w) full(pfn(w));

    function init_collocation()
        tau = full([0, casadi.collocation_points(d, cp)]);
    end

    function init_variables()
        [~, T0, Tf] = ocp.fixed_horizon();
        t0 = yop.cx('t0');
        tf = yop.cx('tf');
        dt = (tf - t0)/N;
        if ocp.n_p > 0
            p  = yop.cx('p');
        else
            p = [];
        end
        init_time();
        init_state();
        init_algebraic();
        init_control();
    end

    function init_time()
        t_n = t0;
        for n=1:N
            t(n) = yop.interpolating_poly(tau, t_n + tau*dt, T0, Tf, N);
            t_n = t_n + dt;
        end
        t(N+1) = yop.interpolating_poly(0, tf, T0, Tf, N);
    end
    
    function init_state()
        for n=1:N
            x_n = yop.cx(['x_' num2str(n)], ocp.n_x, d+1);
            x(n) = yop.interpolating_poly(tau, x_n, T0, Tf, N);
        end
        x_n = yop.cx(['x_' num2str(N+1)], ocp.n_x);
        x(N+1) = yop.interpolating_poly(0, x_n, T0, Tf, N);
    end

    function init_algebraic()
        for n=1:N
            z_n = yop.cx(['z_' num2str(n)], ocp.n_z, d);
            z(n) = yop.interpolating_poly(tau(2:end), z_n, T0, Tf, N);
        end
    end

    function init_control()
        u = yop.interpolating_poly.empty(N, 0);
        for n=1:N
            u_n = yop.cx(['u_' num2str(n)], ocp.n_u);
            u(n) = yop.interpolating_poly(0, u_n, T0, Tf, N);
        end
    end

    function init_derivatives()
        D = yop.interpolating_poly.empty(N, 0);
        for n=1:N
            D(n) = yop.interpolating_poly(tau, [], T0, Tf, N);
        end
    end

    function special_nodes()
        init_derivatives();
        for node = ocp.snodes
            tp_tmp  = [tp; zeros(ocp.n_tp  - length(tp), 1)];
            int_tmp = [I; zeros(ocp.n_int - length(I), 1)];
 
            switch node.type
                case yop.ocp_expr.tp
                    parameterize_timepoint(node, tp_tmp, int_tmp);
                    
                case yop.ocp_expr.int
                    parameterize_integral(node, tp_tmp, int_tmp);
                    
                case yop.ocp_expr.der
                    parameterize_derivative(node, tp_tmp, int_tmp);
            end
        end
    end

    function parameterize_timepoint(ocp_tp, tp_tmp, int_tmp)
        tt = t.value(ocp_tp.timepoint);
        xx = x.value(ocp_tp.timepoint);
        zz = z.value(ocp_tp.timepoint);
        uu = u.value(ocp_tp.timepoint);
        dd = pad(D.value(ocp_tp.timepoint));
        tp = [tp; ocp_tp.fn(t0,tf,tt,xx,zz,uu,p,tp_tmp,int_tmp,dd)];
    end

    function parameterize_integral(i, tp_tmp, int_tmp)
        Ival = 0;
        for n=1:N
            yval = [];
            for r=tau
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = pad(D(n).eval(r));
                yval = [yval, i.fn(...
                    t0, tf, tt, xx, zz, uu, p, tp_tmp, int_tmp, dd)];
            end
            lp = yop.lagrange_polynomial(tau, yval).integrate();
            Ival = Ival + lp.evaluate(1)*dt;
        end
        I = [I; Ival];
    end

    function parameterize_derivative(der, tp_tmp, int_tmp)
        for n=1:N
            yval = [];
            for r=tau
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = pad(D(n).eval(r));
                yval = [yval, der.fn(...
                    t0, tf, tt, xx, zz, uu, p, tp_tmp, int_tmp, dd)];
            end
            lp = yop.lagrange_polynomial(tau, yval);
            yn = lp.differentiate().evaluate(tau)/dt;
            D(n).y = [D(n).y; yn];
        end
    end

    function dd = pad(der)
        dd = [der;  zeros(ocp.n_der  - length(der), 1)];
    end

    function discretize_dynamics()
        for n=1:N % Dynamics
            dx = x(n).differentiate();
            for r=tau(2:end)
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = D(n).eval(r);
                f = ocp.ode.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd);
                a = ocp.alg.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd);
                g = [g; (dx.evaluate(r) - dt*f); a];
            end
        end
        
        for n=1:N % Continuity
            g = [g; x(n).eval(1) - x(n+1).eval(0)];
        end
        
        g_ub = zeros(size(g));
        g_lb = zeros(size(g));
    end

    function pointcons()
        tt = t(1).eval(0);
        xx = x(1).eval(0);
        zz = z(1).eval(0);
        uu = u(1).eval(0);
        dd = D(1).eval(0);
        g = [g; ocp.point.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
        g_ub = [g_ub; ocp.point.ub];
        g_lb = [g_lb; ocp.point.lb];
    end

    function pathcons()
        if ocp.has_path()
            for n=1:N
                tt = t(n).eval(0);
                xx = x(n).eval(0);
                zz = z(n).eval(0);
                uu = u(n).eval(0);
                dd = D(n).eval(0);
                g = [g; ocp.path.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
            end
            tt = t(N+1).eval(0);
            xx = x(N+1).eval(0);
            zz = z(N).eval(1);
            uu = u(N).eval(1);
            dd = D(N).eval(1);
            g = [g; ocp.path.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
            g_ub = [g_ub; repmat(ocp.path.ub, N+1, 1)];
            g_lb = [g_lb; repmat(ocp.path.lb, N+1, 1)];

        end
    end

    function hard_pathcons()
        if ocp.has_hard_path()
            for n=1:N
                for r = tau
                    tt = t(n).eval(r);
                    xx = x(n).eval(r);
                    zz = z(n).eval(r);
                    uu = u(n).eval(r);
                    dd = D(n).eval(r);
                    g = [g; ocp.path_hard.fn(...
                        t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
                end
            end
            tt = t(N+1).eval(0);
            xx = x(N+1).eval(0);
            zz = z(N).eval(1);
            uu = u(N).eval(1);
            dd = D(N).eval(1);
            g = [g; ocp.path_hard.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
            g_ub = [g_ub; repmat(ocp.path_hard.ub, N*(d+1)+1, 1)];
            g_lb = [g_lb; repmat(ocp.path_hard.lb, N*(d+1)+1, 1)];
        end
    end

    function ival_pathcons(expr)
        [I0, If] = ival_bnds(expr);
        
        % Evaluate at the beginning of the interval
        tt = t.value(I0);
        xx = x.value(I0);
        zz = z.value(I0);
        uu = u.value(I0);
        dd = D.value(I0);
        g = [g; expr.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
        g_ub = [g_ub; expr.ub];
        g_lb = [g_lb; expr.lb];
        
        % Evaluate all points within the interval
        [n0, r0, nf, rf] = get_ival_idx(I0, If);
        for n = n0 : min(nf, N) % N+1 is handled by evaluating at If
            for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, d+1)
                tt = t(n).eval(tau(r));
                xx = x(n).eval(tau(r));
                zz = z(n).eval(tau(r));
                uu = u(n).eval(tau(r));
                dd = D(n).eval(tau(r));
                g = [g; expr.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
                g_ub = [g_ub; expr.ub];
                g_lb = [g_lb; expr.lb];
            end
        end
        
        % Evaluate final point
        tt = t.value(If);
        xx = x.value(If);
        zz = z.value(If);
        uu = u.value(If);
        dd = D.value(If);
        g = [g; expr.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd)];
        g_ub = [g_ub; expr.ub];
        g_lb = [g_lb; expr.lb];
    end

    function [I0, If] = ival_bnds(expr)
        I0 = expr.t0;
        If = expr.tf;
        I0 = yop.IF(I0 == yop.initial_timepoint, T0, I0);
        I0 = yop.IF(I0 == yop.final_timepoint  , Tf, I0);
        If = yop.IF(If == yop.initial_timepoint, T0, If);
        If = yop.IF(If == yop.final_timepoint  , Tf, If);
    end

    function [n0, r0, nf, rf] = get_ival_idx(I0, If)
        dT = (Tf-T0)/N;
        
        % First point after I0
        n0 = 1 + floor((I0-T0)/dT);
        if n0 == N+1
            r0 = 1;
        else
            t_n0 = dT*(n0-1);
            tau_I0 = (I0-t_n0)/dT;
            % Does not use >= in order to always pick the next point
            r0 = find(tau - tau_I0 > 0, 1);
            if isempty(r0)
                % Next point is the next control interval
                r0 = 1;
                n0 = n0+1;
            end
        end
        
        % Last point before If
        nf = 1 + floor((If-T0)/dT);
        t_nf = dT*(nf-1);
        tau_If = (If-t_nf)/dT;
        rf = find(tau - tau_If < 0, 1, 'last');
        if isempty(rf)
            if nf==1
                rf = 1;
            else
                rf = length(tau);
                nf = nf-1;
            end
        end
        
    end

    function box_bnd()
        
        dt = (Tf-T0)/N;
        % can be casadi SX/MX, in that case open horizon. If so, it is preferable
        % to use a dummy numeric value as it is faster.
        h = yop.IF(isnumeric(dt), dt, 1);
        
        t0_lb = ocp.t0_lb(T0);
        t0_ub = ocp.t0_ub(T0);
        tf_lb = ocp.tf_lb(T0);
        tf_ub = ocp.tf_ub(T0);
        
        x_ub = ocp.x0_ub(T0);
        x_lb = ocp.x0_lb(T0);
        tt = yop.IF(isnumeric(T0), T0, 1);
        for n=1:N
            for r = tau(yop.IF(n==1,2,1):end)
                x_ub = [x_ub; ocp.x_ub(tt+r)];
                x_lb = [x_lb; ocp.x_lb(tt+r)];
            end
            tt = tt + h;
        end
        x_ub = [x_ub; ocp.xf_ub(Tf)];
        x_lb = [x_lb; ocp.xf_lb(Tf)];
        
        z_ub = [];
        z_lb = [];
        tt = yop.IF(isnumeric(T0), T0, 1);
        for n=1:N
            for r = tau(2:end)
                z_ub = [z_ub; ocp.z_ub(tt+r)];
                z_lb = [z_lb; ocp.z_lb(tt+r)];
            end
            tt = tt + h;
        end
        
        u_ub = ocp.u0_ub(T0);
        u_lb = ocp.u0_lb(T0);
        tt = yop.IF(isnumeric(T0), T0, 1);
        for n=2:N-1
            u_ub = [u_ub; ocp.u_ub(tt)];
            u_lb = [u_lb; ocp.u_lb(tt)];
            tt = tt + h;
        end
        u_ub = [u_ub; ocp.uf_ub(Tf)];
        u_lb = [u_lb; ocp.uf_lb(Tf)];
        
        p_lb = ocp.p_lb(T0);
        p_ub = ocp.p_ub(T0);
        
        w_ub = vertcat(t0_ub, tf_ub, x_ub, z_ub, u_ub, p_ub);
        w_lb = vertcat(t0_lb, tf_lb, x_lb, z_lb, u_lb, p_lb);
    end

end
