classdef ocp_sol < handle
    properties
        independent0
        independentf
        independent
        states
        algebraics
        controls
        parameters
        mx_vars
        t0
        tf
        t
        x
        z
        u
        p
        N
        tau
        dt
    end
    methods
        function obj = ocp_sol(t0, tf, t, x, z, u, p, mx_vars, sol, N, d, cp)
            obj.independent0 = t0;
            obj.independentf = tf;
            obj.independent = t;
            obj.states = x;
            obj.algebraics = z;
            obj.controls = u;
            obj.parameters = p;
            obj.mx_vars = mx_vars;
            obj.t0 = sol.t0;
            obj.tf = sol.tf;
            obj.p = sol.p;
            obj.N = N;
            obj.tau = full([0, casadi.collocation_points(d, cp)]);
            obj.dt = (sol.tf-sol.t0)/N;
            obj.parameterize_polynomials(sol.t, sol.x, sol.z, sol.u);
        end
        
        function n = d(obj)
            % polynomial degree
            n = length(obj.tau) - 1;
        end
        
        function obj = parameterize_polynomials(obj, t, x, z, u)
            ip = @(x,y)yop.interpolating_poly(x, y, obj.t0, obj.tf, obj.N);
            tt = yop.interpolating_poly.empty(obj.N+1, 0);
            xx = yop.interpolating_poly.empty(obj.N+1, 0);
            zz = yop.interpolating_poly.empty(obj.N, 0);
            uu = yop.interpolating_poly.empty(obj.N, 0);
            d = length(obj.tau)-1;
            for n=1:obj.N
                cols = (n-1)*(d+1)+1 : n*(d+1);
                cz   = (n-1)*d+1:n*d;
                tt(n) = ip(obj.tau, t(cols));
                xx(n) = ip(obj.tau, x(:, cols));
                zz(n) = ip(obj.tau(2:end), z(:,cz));
                uu(n) = ip(0, u(:,n));
            end
            tt(n+1) = ip(0, obj.tf);
            xx(n+1) = ip(0, x(:, end));
            obj.t = tt;
            obj.x = xx;
            obj.z = zz;
            obj.u = uu;
        end
        
        function v = value(obj, expr, mag)    
            
            if nargin == 2
                mag = 1;
            end
            
            [vars, tps, ints, ders, sn] = yop.find_special_nodes(expr);
            
            % Ensure that the independent variable always can be plotted by
            % setting its mx value to that of the ocp.
            t_id = 0;
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    % Three is the independent variable. Hard coded, but
                    % efficient.
                    vars{k}.m_value = obj.mx_vars{3};
                    t_id = vars{k}.id;
                end
            end
            
            % Return if the variables thats not know are in the ast
            % This is used in getting values for ocp initial guess, but
            % could also be a user error.
            known_vars = false;
            IDs = [obj.ids, t_id];
            for k=1:length(vars)
                known_vars = known_vars || any(vars{k}.id == IDs);
            end
            if ~known_vars
                v = [];
                return
            end 
            
            % Create functions for the special nodes and expression
            args = { ...
                obj.mx_vars{:}, ...
                tps.mx_vec(), ...
                ints.mx_vec(), ...
                ders.mx_vec() ...
                };
            %             obj.set_mx();
            %             set_mx([tps, ints, ders]);
            for node = [tps, ints, ders]
                mx_expr = fweval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            fn = casadi.Function('fn', args, {fweval(expr)});
            
            % Compute the numerical values of the special nodes
            [tpv, intv, derv] = obj.comp_sn(sn, ...
                n_elem(tps), n_elem(ints), n_elem(ders));
            
            if isa_reducible(expr)
                v = obj.invariant_value(fn, tpv, intv, derv);
            elseif is_ival(expr)
                v = obj.interval_value(expr, fn, tpv, intv, derv, mag);
            else
                v = obj.variant_value(fn, tpv, intv, derv, mag);
            end
        end
        
        function IDs = ids(obj)
            IDs = [obj.independent0.ids, obj.independentf.ids, ...
                obj.independent.ids, obj.states.ids, obj.algebraics.ids, ...
                obj.controls.ids, obj.parameters.ids];
        end
        
        function v = invariant_value(obj, expr, tps, ints, ders)
            v = full(expr(obj.t0, obj.tf, ...
                obj.t(1).evaluate(0), ...
                obj.x(1).evaluate(0), ...
                obj.z(1).evaluate(0), ...
                obj.u(1).evaluate(0), ...
                obj.p, tps, ints, ders(1).evaluate(0)));
        end
        
        function v = variant_value(obj, expr, tps, ints, ders, mag)
            tt=[]; xx=[]; zz=[]; uu=[]; dd=[];
            for n=1:obj.N
                un = obj.u(n).y;
                for r=1:length(obj.tau)-1
                    tt = [tt, obj.t(n).y(r)];
                    xx = [xx, obj.x(n).y(:,r)];
                    zz = [zz, obj.z(n).evaluate(obj.tau(r))];
                    uu = [uu, un];
                    dd = [dd, ders(n).evaluate(obj.tau(r))];
                    dT = obj.tau(r+1)-obj.tau(r);
                    for k=1:mag-1 % Magnification
                        tau_k = obj.tau(r) + k/mag*dT;
                        tt = [tt, obj.t(n).evaluate(tau_k)];
                        xx = [xx, obj.x(n).evaluate(tau_k)];
                        zz = [zz, obj.z(n).evaluate(tau_k)];
                        uu = [uu, un];
                        dd = [dd, ders(n).evaluate(tau_k)];
                    end
                end
                tt = [tt, obj.t(n).y(r+1)];
                xx = [xx, obj.x(n).y(:,r+1)];
                zz = [zz, obj.z(n).y(:,r)];
                uu = [uu, un];
                dd = [dd, ders(n).evaluate(obj.tau(r+1))];
            end
            tt = [tt, obj.t(n+1).y];
            xx = [xx, obj.x(n+1).y(:)];
            zz = [zz, obj.z(n).evaluate(1)];
            uu = [uu, obj.u(n).y];
            dd = [dd, ders(n).evaluate(1)];
            v = full(expr(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
        end
        
        function v = interval_value(obj, expr, fn, tps, ints, ders, mag)
            [I0, If] = get_ival(expr);
            I0 = yop.IF(I0==yop.initial_timepoint, obj.t0, I0);
            I0 = yop.IF(I0==yop.final_timepoint  , obj.tf, I0);
            If = yop.IF(If==yop.initial_timepoint, obj.t0, If);
            If = yop.IF(If==yop.final_timepoint  , obj.tf, If);
            
            % Evaluate at the beginning of the interval
            tt = obj.t.value(I0); 
            xx = obj.x.value(I0);
            zz = obj.z.value(I0);
            uu = obj.u.value(I0);
            dd = ders.value(I0); 
            
            % Evaluate all points within the interval
            [n0, r0, nf, rf] = obj.get_ival_idx(I0, If);
            for n=n0:min(nf, obj.N)
                un = obj.u(n).y;
                for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, obj.d)
                    tt = [tt, obj.t(n).y(r)];
                    xx = [xx, obj.x(n).y(:, r)];
                    zz = [zz, obj.z(n).evaluate(obj.tau(r))];
                    uu = [uu, un];
                    dd = [dd, ders(n).evaluate(obj.tau(r))];
                    if n==nf && r==rf
                        t_nf = obj.dt*(nf-1);
                        tau_If = (If-t_nf)/obj.dt;
                        dT = tau_If - obj.tau(r);
                    else
                        dT = obj.tau(r+1)-obj.tau(r);
                    end
                    for k=1:mag-1 % Magnification
                        tau_k = obj.tau(r) + k/mag*dT;
                        tt = [tt, obj.t(n).evaluate(tau_k)];
                        xx = [xx, obj.x(n).evaluate(tau_k)];
                        zz = [zz, obj.z(n).evaluate(tau_k)];
                        uu = [uu, un];
                        dd = [dd, ders(n).evaluate(tau_k)];
                    end
                end
                % Last collocation point is evaluated here
                if n ~= nf
                    tt = [tt, obj.t(n).y(r+1)];
                    xx = [xx, obj.x(n).y(:,r+1)];
                    zz = [zz, obj.z(n).y(:,r)];
                    uu = [uu, obj.u(n).evaluate(obj.tau(r))];
                    dd = [dd, ders(n).evaluate(obj.tau(r+1))];
                end
            end
            
            % Evaluate last point of interval
            tt = [tt, obj.t.value(If)]; 
            xx = [xx, obj.x.value(If)];
            zz = [zz, obj.z.value(If)];
            uu = [uu, obj.u.value(If)];
            dd = [dd, ders.value(If)]; 
            
            v = full(fn(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
        end
        
        function [n0, r0, nf, rf] = get_ival_idx(obj, I0, If)
            
            % First point after I0
            n0 = 1 + floor((I0-obj.t0)/obj.dt);
            if n0 == obj.N+1
                r0 = 1;
            else
                t_n0 = obj.dt*(n0-1);
                tau_I0 = (I0-t_n0)/obj.dt;
                % Does not use >= in order to always pick the next point
                r0 = find(obj.tau - tau_I0 > 0, 1);
                if isempty(r0)
                    r0 = 1;
                    n0 = n0+1;
                end
            end
            
            % Last point before If
            nf = 1 + floor((If-obj.t0)/obj.dt);
            t_nf = obj.dt*(nf-1);
            tau_If = (If-t_nf)/obj.dt;
            rf = find(obj.tau - tau_If < 0, 1, 'last');
            if isempty(rf)
                if nf==1
                    rf = 1;
                else
                    rf = length(obj.tau);
                    nf = nf-1;
                end
            end
            
        end
        
        function [tps, ints, ders] = comp_sn(obj, sn, n_tp, n_int, n_der)
            tps = [];
            ints = [];
            ders = yop.interpolating_poly.empty(obj.N, 0);
            for n=1:obj.N
                ders(n) = yop.interpolating_poly(obj.tau, [], obj.t0, ...
                    obj.tf, obj.N);
            end
            for node = sn
                tmp_tp  = [tps;  zeros(n_tp  - length(tps), 1)];
                tmp_int = [ints; zeros(n_int - length(ints), 1)];
                switch node.type
                    case yop.ocp_expr.tp
                        tp = obj.compute_timepoint( ...
                            node, tmp_tp, tmp_int, ders, n_der);
                        tps = [tps; tp(:)];
                        
                    case yop.ocp_expr.int
                        int = obj.compute_integral( ...
                            node, tmp_tp, tmp_int, ders, n_der);
                        ints = [ints; int(:)];
                        
                    case yop.ocp_expr.der
                        ders = obj.compute_derivative(node, tmp_tp, ...
                            tmp_int, ders, n_der);
                        
                    otherwise
                        error(yop.msg.unexpected_error);
                end
            end
        end
        
        function val = compute_timepoint(obj, tp, tps, ints, ders, n_der)
            dd = ders.value(tp.timepoint);
            dd = [dd;  zeros(n_der  - length(dd), 1)];
            val = full(tp.fn(obj.t0, obj.tf, ...
                obj.t.value(tp.timepoint), ...
                obj.x.value(tp.timepoint), ...
                obj.z.value(tp.timepoint), ...
                obj.u.value(tp.timepoint), ...
                obj.p, tps, ints, dd));
        end
        
        function I = compute_integral(obj, int, tps, ints, ders, n_der)
            I = 0;
            for n=1:obj.N
                yval = [];
                for r=1:length(obj.tau)
                    dd = ders(n).evaluate(obj.tau(r));
                    dd = [dd;  zeros(n_der  - length(dd), 1)];
                    val_r = full(int.fn(obj.t0, obj.tf, ...
                        obj.t(n).y(r), ...
                        obj.x(n).y(:, r), ...
                        obj.z(n).evaluate(obj.tau(r)), ...
                        obj.u(n).y, ...
                        obj.p, tps, ints, dd));
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(obj.tau, yval).integrate();
                I = I + lp.evaluate(1)*obj.dt;
            end
        end
        
        function ders = compute_derivative(obj,der,tps,ints,ders,n_der)
            for n=1:obj.N
                yval = [];
                for r=1:length(obj.tau)
                    tt = obj.t(n).y(r);
                    xx = obj.x(n).y(:, r);
                    zz = obj.z(n).evaluate(obj.tau(r));
                    uu = obj.u(n).y;
                    pp = obj.p;
                    dd = ders(n).evaluate(obj.tau(r));
                    dd = [dd;  zeros(n_der  - length(dd), 1)];
                    val_r = full(der.fn(obj.t0, obj.tf, tt, xx, zz, uu, pp, ...
                        tps, ints, dd));
                    yval = [yval, val_r(:)];
                end
                yn = yop.lagrange_polynomial(obj.tau, ...
                    yval).differentiate().evaluate(obj.tau)/obj.dt;
                ders(n).y = [ders(n).y; yn];
            end
        end
        
        function args = filter(obj, input)
            args = {};
            
            % First filter out magnification
            mag = 1;
            k = 1;
            while k <= length(input)
                if strcmp(input{k}, 'mag')
                    mag = input{k+1};
                    k = k + 2;
                else
                    args{end+1} = input{k};
                    k = k + 1;
                end
            end
            
            % Get values
            for k=1:length(args)
                if isa(args{k}, 'yop.ast_node')
                    args{k} = obj.value(args{k}, round(mag));
                end
            end
        end
        
        %         function set_mx(obj)
        %             obj.independent0.set_mx();
        %             obj.independentf.set_mx();
        %             obj.independent.set_mx();
        %             obj.states.set_mx();
        %             obj.algebraics.set_mx();
        %             obj.controls.set_mx();
        %             obj.parameters.set_mx();
        %         end
        
        
        function varargout = plot(obj, varargin)
            args = obj.filter(varargin);
            h = plot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = plot3(obj, varargin)
            args = obj.filter(varargin);
            h = plot3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = stairs(obj, varargin)
            args = obj.filter(varargin);
            h = stairs(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = errorbar(obj, varargin)
            args = obj.filter(varargin);
            h = errorbar(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = area(obj, varargin)
            args = obj.filter(varargin);
            h = area(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = stackedplot(obj, varargin)
            args = obj.filter(varargin);
            h = stackedplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = loglog(obj, varargin)
            args = obj.filter(varargin);
            h = loglog(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = semilogx(obj, varargin)
            args = obj.filter(varargin);
            h = semilogx(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = semilogy(obj, varargin)
            args = obj.filter(varargin);
            h = semilogy(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fplot(obj, varargin)
            args = obj.filter(varargin);
            h = fplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fplot3(obj, varargin)
            args = obj.filter(varargin);
            h = fplot3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fimplicit(obj, varargin)
            args = obj.filter(varargin);
            h = fimplicit(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = scatter(obj, varargin)
            args = obj.filter(varargin);
            h = scatter(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = scatter3(obj, varargin)
            args = obj.filter(varargin);
            h = scatter3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bubblechart(obj, varargin)
            args = obj.filter(varargin);
            h = bubblechart(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bubblechart3(obj, varargin)
            args = obj.filter(varargin);
            h = bubblechart3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = swarmchart(obj, varargin)
            args = obj.filter(varargin);
            h = swarmchart(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = swarmchart3(obj, varargin)
            args = obj.filter(varargin);
            h = swarmchart3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = spy(obj, varargin)
            args = obj.filter(varargin);
            h = spy(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = histogram(obj, varargin)
            args = obj.filter(varargin);
            h = histogram(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = histogram2(obj, varargin)
            args = obj.filter(varargin);
            h = histogram2(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = pie(obj, varargin)
            args = obj.filter(varargin);
            h = pie(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        
        function varargout = pie3(obj, varargin)
            args = obj.filter(varargin);
            h = pie3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = scatterhistogram(obj, varargin)
            args = obj.filter(varargin);
            h = scatterhistogram(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = wordcloud(obj, varargin)
            args = obj.filter(varargin);
            h = wordcloud(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bubblecloud(obj, varargin)
            args = obj.filter(varargin);
            h = bubblecloud(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = heatmap(obj, varargin)
            args = obj.filter(varargin);
            h = heatmap(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = parllelplot(obj, varargin)
            args = obj.filter(varargin);
            h = parllelplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = plotmatrix(obj, varargin)
            args = obj.filter(varargin);
            h = plotmatrix(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bar(obj, varargin)
            args = obj.filter(varargin);
            h = bar(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = barh(obj, varargin)
            args = obj.filter(varargin);
            h = barh(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bar3(obj, varargin)
            args = obj.filter(varargin);
            h = bar3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = bar3h(obj, varargin)
            args = obj.filter(varargin);
            h = bar3h(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = pareto(obj, varargin)
            args = obj.filter(varargin);
            h = pareto(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = stem(obj, varargin)
            args = obj.filter(varargin);
            h = stem(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = stem3(obj, varargin)
            args = obj.filter(varargin);
            h = stem3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = geoplot(obj, varargin)
            args = obj.filter(varargin);
            h = geoplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = geoscatter(obj, varargin)
            args = obj.filter(varargin);
            h = geoscatter(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = geobubble(obj, varargin)
            args = obj.filter(varargin);
            h = geobubble(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = polarplot(obj, varargin)
            args = obj.filter(varargin);
            h = polarplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = polarhistogram(obj, varargin)
            args = obj.filter(varargin);
            h = polarhistogram(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = polarscatter(obj, varargin)
            args = obj.filter(varargin);
            h = polarscatter(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = polarbubblechart(obj, varargin)
            args = obj.filter(varargin);
            h = polarbubblechart(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = compass(obj, varargin)
            args = obj.filter(varargin);
            h = compass(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = ezpolar(obj, varargin)
            args = obj.filter(varargin);
            h = ezpolar(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = contour(obj, varargin)
            args = obj.filter(varargin);
            h = contour(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = contourf(obj, varargin)
            args = obj.filter(varargin);
            h = contourf(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = contour3(obj, varargin)
            args = obj.filter(varargin);
            h = contour3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = contourslice(obj, varargin)
            args = obj.filter(varargin);
            h = contourslice(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fcontour(obj, varargin)
            args = obj.filter(varargin);
            h = fcontour(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = quiver(obj, varargin)
            args = obj.filter(varargin);
            h = quiver(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = quiver3(obj, varargin)
            args = obj.filter(varargin);
            h = quiver3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = feather(obj, varargin)
            args = obj.filter(varargin);
            h = feather(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = surf(obj, varargin)
            args = obj.filter(varargin);
            h = surf(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = surfc(obj, varargin)
            args = obj.filter(varargin);
            h = surfc(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = surfl(obj, varargin)
            args = obj.filter(varargin);
            h = surfl(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = ribbon(obj, varargin)
            args = obj.filter(varargin);
            h = ribbon(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = pcolor(obj, varargin)
            args = obj.filter(varargin);
            h = surf(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fsurf(obj, varargin)
            args = obj.filter(varargin);
            h = fsurf(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fimplicit3(obj, varargin)
            args = obj.filter(varargin);
            h = fimplicit3(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = mesh(obj, varargin)
            args = obj.filter(varargin);
            h = mesh(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = meshc(obj, varargin)
            args = obj.filter(varargin);
            h = meshc(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = meshz(obj, varargin)
            args = obj.filter(varargin);
            h = meshz(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = waterfall(obj, varargin)
            args = obj.filter(varargin);
            h = waterfall(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = fmesh(obj, varargin)
            args = obj.filter(varargin);
            h = fmesh(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = streamline(obj, varargin)
            args = obj.filter(varargin);
            h = streamline(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = streamslice(obj, varargin)
            args = obj.filter(varargin);
            h = streamslice(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = streamparticles(obj, varargin)
            args = obj.filter(varargin);
            h = streamparticles(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = streamribbon(obj, varargin)
            args = obj.filter(varargin);
            h = streamribbon(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = streamtube(obj, varargin)
            args = obj.filter(varargin);
            h = streamtube(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = coneplot(obj, varargin)
            args = obj.filter(varargin);
            h = coneplot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = slice(obj, varargin)
            args = obj.filter(varargin);
            h = slice(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = image(obj, varargin)
            args = obj.filter(varargin);
            h = image(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = imagesc(obj, varargin)
            args = obj.filter(varargin);
            h = imagesc(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
    end
end