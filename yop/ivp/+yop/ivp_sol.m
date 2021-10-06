classdef ivp_sol < handle
    
    properties
        sol
        solver
        ivp
        t0
        tf
        t
        x
        z
        p
        N
        d
        tau
    end
    
    methods
        
        function obj = ivp_sol(ivp, sol, solver, t, p, N, tau)
            obj.sol = sol;
            obj.solver = solver;
            obj.ivp = ivp;
            obj.t0 = t(1);
            obj.tf = t(end);
            obj.N = N;
            obj.tau = tau;
            
            d = length(tau)-1;
            x = full(sol.xf);
            z = full(sol.zf);
            p = p;
            ts = {};
            xs = {};
            zs = {};
            for n=1:N
                cols = (n-1)*(d+1)+1 : n*(d+1);
                ts{n} = t(cols);
                xs{n} = x(:, cols);
                zs{n} = z(:, cols);
            end
            ts{N+1} = t(end);
            xs{N+1} = x(:,end);
            zs{N+1} = z(:,end);
            
            tlp = yop.sim_poly.empty(1,0);
            xlp = yop.sim_poly.empty(1,0);
            zlp = yop.sim_poly.empty(1,0);
            for n=1:N
                T0 = ts{n}(1);
                Tf = ts{n+1}(1);
                dt = (Tf-T0)/N;
                tlp(n) = yop.sim_poly(tau, ts{n}, T0, Tf);
                xlp(n) = yop.sim_poly(tau, xs{n}, T0, Tf);
                zlp(n) = yop.sim_poly(tau, zs{n}, T0, Tf);
            end
            tlp(N+1) = yop.sim_poly(0, ts{N+1}, Tf, Tf);
            xlp(N+1) = yop.sim_poly(0, xs{N+1}, Tf, Tf);
            zlp(N+1) = yop.sim_poly(0, zs{N+1}, Tf, Tf);
            
            obj.t = tlp;
            obj.x = xlp;
            obj.z = zlp;
            obj.p = p;
        end
        
        function v = value(obj, expr)
            expr = yop.ocp_expr(expr);
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr.ast);
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.ivp.independent.mx;
                end
            end
            
            args = { ...
                mx_vec(obj.ivp.independent), ...
                mx_vec(obj.ivp.states), ...
                mx_vec(obj.ivp.algebraics), ...
                mx_vec(obj.ivp.parameters) ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx([obj.ivp.independent, obj.ivp.states, ...
                obj.ivp.algebraics, obj.ivp.parameters]);
            set_mx([tps, ints, ders]);
            
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});
            
            [tpv, intv] = obj.parameterize_special_nodes(sn, n_elem(tps), n_elem(ints));
            
            if expr.is_transcription_invariant
                v = expr.fn( ...
                    obj.t(1).evaluate(0), ...
                    obj.x(1).evaluate(0), ...
                    obj.z(1).evaluate(0), ...
                    obj.p, tpv, intv, []);
            else
                tt = mat(obj.t);
                xx = mat(obj.x);
                zz = mat(obj.z);
                v = expr.fn(tt,xx,zz,obj.p,tpv,intv,[]);
            end
            v = full(v);
        end
        
        function [tpv, intv] = parameterize_special_nodes(obj, sn, n_tp, n_int)
            tpv = [];
            intv = [];
            ders = [];
            for node = sn
                tmp_tp  = [tpv;  zeros(n_tp  - length(tpv), 1)];
                tmp_int = [intv; zeros(n_int - length(intv), 1)];
                switch node.type
                    case yop.ocp_expr.tp
                        tp = obj.compute_timepoint( ...
                            node, tmp_tp, tmp_int, ders);
                        tpv = [tpv; tp(:)];
                        
                    case yop.ocp_expr.int
                        int = obj.compute_integral( ...
                            node, tmp_tp, tmp_int, ders);
                        intv = [intv; int(:)];
                        
                    case yop.ocp_expr.der
                        error(yop.msg.not_implemented);
                        warning(['When implementing do not forget to ' ...
                            'multiply derivative with step length']);
                    otherwise
                        error(yop.msg.unexpected_error);
                end
            end
        end
        
        function v = compute_timepoint(obj, tp, tpv, intv, der)
            v = tp.fn( ...
                obj.t.value(tp.timepoint), ...
                obj.x.value(tp.timepoint), ...
                obj.z.value(tp.timepoint), ...
                obj.p, tpv, intv, der);
        end
        
        function I = compute_integral(obj, int, tpv, intv, der)
            I = 0;
            for n=1:obj.N
                yval = [];
                for r=obj.tau
                    tt = obj.t(n).evaluate(r);
                    xx = obj.x(n).evaluate(r);
                    zz = obj.z(n).evaluate(r);
                    val_r = int.fn(tt, xx, zz, obj.p, tpv, intv, der);
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(obj.tau, yval).integrate();
                dt = obj.t(n).tf - obj.t(n).t0;
                I = I + lp.evaluate(1)*dt;
            end
        end
        
        
        function args = filter(obj, input)           
            for k=1:length(input)
                if isa(input{k}, 'yop.node')
                    args{k} = obj.value(input{k});
                else
                    args{k} = input{k};
                end
            end
        end
        
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