function nlp = direct_collocation(ocp, N, dx, cpx, du, cpu)

                       yop.progress.nlp_building();
T0=[]; Tf=[]; t=[]; t0=[]; tf=[]; p=[]; dt=[]; taux=[]; tp=[]; I=[]; D=[]; 
g=[]; g_ub=[]; g_lb=[]; w_ub=[]; w_lb=[]; w0=[]; h=[]; tauu=[];
t = yop.interpolating_poly.empty(0, N+1);
x = yop.interpolating_poly.empty(0, N+1);
z = yop.interpolating_poly.empty(0, N);
u = yop.interpolating_poly.empty(0, N);
D = yop.interpolating_poly.empty(0, N);

init_collocation();
init_variables();      yop.progress.nlp_initialized();
special_nodes();       yop.progress.nlp_special_nodes_done();
discretize_dynamics(); yop.progress.nlp_dynamics_done();
pointcons();           yop.progress.nlp_pointcons_done();
pathcons();            
hard_pathcons();
for pk = ocp.path_ival
    ival_pathcons(pk);
end
continuous_control();
                       yop.progress.nlp_pathcons_done();

w = vertcat(t0, tf, x.vec(), z.vec(), u.vec(), p);
box_bnd();             yop.progress.nlp_box_done();
initial_guess();       yop.progress.nlp_guess_done();

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
nlp.w0 = w0;
nlp.g = g;
nlp.g_ub = g_ub;
nlp.g_lb = g_lb;
nlp.t0 = t0;
nlp.tf = tf;
nlp.t  = t;
nlp.x  = x;
nlp.z  = z;
nlp.u  = u;
nlp.p  = p;
nlp.ocp_t0 = t0fn;
nlp.ocp_tf = tffn;
nlp.ocp_t  = tfn;
nlp.ocp_x  = xfn;
nlp.ocp_z  = zfn;
nlp.ocp_u  = ufn;
nlp.ocp_p  = pfn;  yop.progress.nlp_completed();

    function init_collocation()
        taux = full([0, casadi.collocation_points(dx, cpx)]);
        tauu = full(casadi.collocation_points(du, cpu));
    end

    function init_variables()
        [T0, Tf] = ocp.get_horizon();
        h = (Tf-T0)/N;
        t0 = yop.cx('t0');
        tf = yop.cx('tf');
        dt = (ocp.descale_tf(tf) - ocp.descale_t0(t0))/N;
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
        T = T0;
        for n=1:N
            t(n) = yop.interpolating_poly(taux, t_n + taux*dt, T, T+h);
            t_n = t_n + dt;
            T = T + h;
        end
        assert(yop.EQ(T, Tf));
        t(N+1) = yop.interpolating_poly(0, tf, Tf, Tf);
    end
    
    function init_state()
        T = T0;
        for n=1:N
            x_n = yop.cx(['x_' num2str(n)], ocp.n_x, dx+1);
            x(n) = yop.interpolating_poly(taux, x_n, T, T+h);
            T = T + h;
        end
        assert(yop.EQ(T, Tf));
        x_n = yop.cx(['x_' num2str(N+1)], ocp.n_x);
        x(N+1) = yop.interpolating_poly(0, x_n, Tf, Tf);
    end

    function init_algebraic()
        T = T0;
        for n=1:N
            z_n = yop.cx(['z_' num2str(n)], ocp.n_z, dx);
            z(n) = yop.interpolating_poly(taux(2:end), z_n, T, T+h);
            T = T + h;
        end
        assert(yop.EQ(T, Tf));
    end

    function init_control()
        T = T0;
        for n=1:N
            u_n = yop.cx(['u_' num2str(n)], ocp.n_u, du);
            u(n) = yop.interpolating_poly(tauu, u_n, T, T+h);
            T = T + h;
        end
        assert(yop.EQ(T, Tf));
    end

    function continuous_control()
        for n=1:N-1
            c = u(n).eval(1)-u(n+1).eval(0);
            g = [g; c];
            g_ub = [g_ub; zeros(size(c))];
            g_lb = [g_lb; zeros(size(c))];
        end
    end

    function init_derivatives()
        T = T0;
        for n=1:N
            D(n) = yop.interpolating_poly(taux, [], T, T+h);
            T = T + h;
        end
        assert(yop.EQ(T, Tf));
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
            for r=taux
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = pad(D(n).eval(r));
                yval = [yval, i.fn(...
                    t0, tf, tt, xx, zz, uu, p, tp_tmp, int_tmp, dd)];
            end
            lp = yop.lagrange_polynomial(taux, yval).integrate();
            Ival = Ival + lp.evaluate(1)*dt;
        end
        I = [I; Ival];
    end

    function parameterize_derivative(der, tp_tmp, int_tmp)
        for n=1:N
            yval = [];
            for r=taux
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = pad(D(n).eval(r));
                yval = [yval, der.fn(...
                    t0, tf, tt, xx, zz, uu, p, tp_tmp, int_tmp, dd)];
            end
            lp = yop.lagrange_polynomial(taux, yval);
            yn = lp.differentiate().evaluate(taux)/dt;
            D(n).y = [D(n).y; yn];
        end
    end

    function dd = pad(der)
        dd = [der;  zeros(ocp.n_der  - length(der), 1)];
    end

    function discretize_dynamics()
        for n=1:N % Dynamics
            derx = x(n).differentiate();
            for r=taux(2:end)
                tt = t(n).eval(r);
                xx = x(n).eval(r);
                zz = z(n).eval(r);
                uu = u(n).eval(r);
                dd = D(n).eval(r);
                f = ocp.ode.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd);
                a = ocp.alg.fn(t0, tf, tt, xx, zz, uu, p, tp, I, dd);
                g = [g; (derx.evaluate(r) - dt*f); a];
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
                for r = taux
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
            g_ub = [g_ub; repmat(ocp.path_hard.ub, N*(dx+1)+1, 1)];
            g_lb = [g_lb; repmat(ocp.path_hard.lb, N*(dx+1)+1, 1)];
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
            for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, dx+1)
                tt = t(n).eval(taux(r));
                xx = x(n).eval(taux(r));
                zz = z(n).eval(taux(r));
                uu = u(n).eval(taux(r));
                dd = D(n).eval(taux(r));
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
        % First point after I0
        n0 = 1 + floor((I0-T0)/h);
        if n0 == N+1
            r0 = 1;
        else
            t_n0 = h*(n0-1);
            tau_I0 = (I0-t_n0)/h;
            % Does not use >= in order to always pick the next point
            r0 = find(taux - tau_I0 > 0, 1);
            if isempty(r0)
                % Next point is the next control interval
                r0 = 1;
                n0 = n0+1;
            end
        end
        
        % Last point before If
        nf = 1 + floor((If-T0)/h);
        t_nf = h*(nf-1);
        tau_If = (If-t_nf)/h;
        rf = find(taux - tau_If < 0, 1, 'last');
        if isempty(rf)
            if nf==1
                rf = 1;
            else
                rf = length(taux);
                nf = nf-1;
            end
        end
        
    end

    function box_bnd()
        t0_lb = ocp.t0_lb(T0);
        t0_ub = ocp.t0_ub(T0);
        tf_lb = ocp.tf_lb(T0);
        tf_ub = ocp.tf_ub(T0);
        
        x_ub = ocp.x0_ub(T0);
        x_lb = ocp.x0_lb(T0);
        T = T0;
        for n=1:N
            for r = taux(yop.IF(n==1,2,1):end)*h
                x_ub = [x_ub; ocp.x_ub(T+r)];
                x_lb = [x_lb; ocp.x_lb(T+r)];
            end
            T = T + h;
        end
        x_ub = [x_ub; ocp.xf_ub(Tf)];
        x_lb = [x_lb; ocp.xf_lb(Tf)];
        
        z_ub = [];
        z_lb = [];
        T = T0;
        for n=1:N
            for r = taux(2:end)*h
                z_ub = [z_ub; ocp.z_ub(T+r)];
                z_lb = [z_lb; ocp.z_lb(T+r)];
            end
            T = T + h;
        end
        
        
        % First interval
        if du == 1
            u_ub = ocp.u0_ub(T0);
            u_lb = ocp.u0_lb(T0);
        else
            g    = [g   ; u(1).eval(0)];
            g_ub = [g_ub; ocp.u0_ub(T0)];
            g_lb = [g_lb; ocp.u0_lb(T0)];
            
            u_ub=[]; u_lb=[];
            for r = tauu*h
                u_ub = [u_ub; ocp.u_ub(r)];
                u_lb = [u_lb; ocp.u_lb(r)];
            end
        end
        
        % Middle intervals
        T = T0 + h;
        for n=2:N-1
            for r = tauu*h
                u_ub = [u_ub; ocp.u_ub(T+r)];
                u_lb = [u_lb; ocp.u_lb(T+r)];
            end
            T = T + h;
        end
        
        % Last interval
        if du == 1    
            u_ub = [u_ub; ocp.uf_ub(Tf)];
            u_lb = [u_lb; ocp.uf_lb(Tf)];
        else
            for r = tauu(1:end-1)*h
                u_ub = [u_ub; ocp.u_ub(T+r)];
                u_lb = [u_lb; ocp.u_lb(T+r)];
            end
            if ~strcmp(cpu, 'radau')
                % Final point is not a nlp variable, so we add final box
                % constraint as a nonlinear constraint.
                
                % First evaluate last collocation point
                % This is done here, since radau points uses uf_ub instead
                % of u_ub for last point.
                u_ub = [u_ub; ocp.u_ub(T+tauu(end)*h)];
                u_lb = [u_lb; ocp.u_lb(T+tauu(end)*h)];
                
                g    = [g   ; u(end).eval(1)];
                g_ub = [g_ub; ocp.uf_ub(Tf)];
                g_lb = [g_lb; ocp.uf_lb(Tf)];
            else
                u_ub = [u_ub; ocp.uf_ub(Tf)];
                u_lb = [u_lb; ocp.uf_lb(Tf)];
                
            end
        end
        
        %         u_ub = ocp.u0_ub(T0);
        %         u_lb = ocp.u0_lb(T0);
        %         T = T0;
        %         for n=2:N-1
        %             u_ub = [u_ub; ocp.u_ub(T)];
        %             u_lb = [u_lb; ocp.u_lb(T)];
        %             T = T + h;
        %         end
        %         u_ub = [u_ub; ocp.uf_ub(Tf)];
        %         u_lb = [u_lb; ocp.uf_lb(Tf)];
        
        p_lb = ocp.p_lb(T0);
        p_ub = ocp.p_ub(T0);
        
        w_ub = vertcat(t0_ub, tf_ub, x_ub, z_ub, u_ub, p_ub);
        w_lb = vertcat(t0_lb, tf_lb, x_lb, z_lb, u_lb, p_lb);
    end

    function initial_guess()
        if ocp.has_initial_guess()
            [t0_0, tf_0, t_0, x_0, z_0, u_0, p_0] = ocp.initial_guess();
            [t_x, t_z, t_u] = grid_at_iv(t_0(1), t_0(end));
            x0 = interpolate(t_0, x_0, t_x(:))';
            z0 = interpolate(t_0, z_0, t_z(:))';
            u0 = interpolate(t_0, u_0, t_u(:))';
            w0 = [t0_0; tf_0; x0(:); z0(:); u0(:); p_0];
        else
            w0 = ones(size(w));
        end
    end
    
    function [tx, tz, tu] = grid_at_iv(t0_0, tf_0)
        % Time points for the initial values
        tx=[]; tz=[]; tu=[];
        h0 = (tf_0-t0_0)/N;
        for k=1:N
            tx = [tx, t0_0 + taux*h0];
            tz = [tz, t0_0 + taux(2:end)*h0];
            tu = [tu, t0_0 + tauu*h0];
            %tu = [tu, t0_0];
            t0_0 = t0_0 + h0;
        end
        tx(end+1) = tf_0;
    end

end

function v = interpolate(x_data, y_data, x_interp)
if isempty(y_data)
    v = [];
    return;
end
v = interp1(x_data, y_data, x_interp);
end
