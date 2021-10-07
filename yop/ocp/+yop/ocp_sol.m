classdef ocp_sol < handle
    properties
        ocp_vars
        mx_args
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
        function obj = ocp_sol(ocp_vars,mx_args,t0,tf,t,x,z,u,p,N,d,cp)
            obj.ocp_vars = ocp_vars;
            obj.mx_args = mx_args;
            obj.t0 = t0;
            obj.tf = tf;
            obj.p = p;
            obj.N = N;
            obj.tau = full([0, casadi.collocation_points(d, cp)]);
            obj.dt = (tf-t0)/N;
            obj.parameterize_polynomials(t, x, z, u);
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
            tt(obj.N+1) = ip(0, obj.tf);
            xx(obj.N+1) = ip(0, x(:, end));
            obj.t = tt;
            obj.x = xx;
            obj.z = zz;
            obj.u = uu;
        end
        
        function v = value(obj, expr, mag)            
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr);
            
            % Ensure that the independent variable always can be plotted by
            % setting its mx value to that of the ocp.
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    % Three is the independent variable. Hard coded but
                    % efficient.
                    vars{k}.m_value = obj.ocp_vars(3).mx;
                end
            end
            
            % Create functions for the special nodes and expression
            args = {obj.mx_args{:}, mx_vec(tps), mx_vec(ints), ...
                    mx_vec(ders)};
            set_mx(obj.ocp_vars);
            set_mx([tps, ints, ders]);
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            fn = casadi.Function('fn', args, {fw_eval(expr)});
            
            % Compute the numerical values of the special nodes
            [tpv, intv, derv] = obj.comp_sn(sn, ...
                n_elem(tps), n_elem(ints), n_elem(ders));
            
            if is_transcription_invariant(expr)
                v = fn(obj.t0, obj.tf, ...
                    obj.t(1).evaluate(0), ...
                    obj.x(1).evaluate(0), ...
                    obj.z(1).evaluate(0), ...
                    obj.u(1).evaluate(0), ...
                    obj.p, tpv, intv, derv);
            elseif mag > 1
                teval = [];
                tvec = vec(obj.t);
                for n=1:length(tvec)-1
                    teval(end+1) = tvec(n);
                    T = tvec(n+1) - tvec(n);
                    for k=1:mag-1
                        teval(end+1) = tvec(n) + T*k/mag;
                    end
                end
                teval(end+1) = tvec(end);
                
                v = [];
                for tp=teval
                    vk = fn(obj.t0, obj.tf, ...
                        obj.t.value(tp), ...
                        obj.x.value(tp), ...
                        obj.z.value(tp), ...
                        obj.u.value(tp), ...
                        obj.p, tpv, intv, derv);
                    v = [v, vk];
                end
                
            else % eval at all collocation points
                v = [];
                for n=1:obj.N
                    uu = obj.u(n).y;
                    v = [v, fn(obj.t0, obj.tf, ...
                        obj.t(n).y(1), ...
                        obj.x(n).y(:,1), ...
                        obj.z(n).evaluate(obj.tau(1)), ...
                        uu, ...
                        obj.p, tpv, intv, derv)];
                    for r=2:length(obj.tau)
                        v = [v, fn(obj.t0, obj.tf, ...
                        obj.t(n).y(r), ...
                        obj.x(n).y(:,r), ...
                        obj.z(n).y(:,r-1), ...
                        uu, ...
                        obj.p, tpv, intv, derv)];
                    end
                end
                v = [v, fn(obj.t0, obj.tf, ...
                    obj.t(n+1).y(1), ...
                    obj.x(n+1).y(:), ...
                    obj.z(n).evaluate(obj.tau(1)), ...
                    uu, ...
                    obj.p, tpv, intv, derv)];
            end
            v = full(v);
        end
        
        function [tps, ints, ders] = comp_sn(obj, sn, n_tp, n_int, n_der)
            tps = [];
            ints = [];
            ders = [];
            for node = sn
                tmp_tp  = [tps;  zeros(n_tp  - length(tps), 1)];
                tmp_int = [ints; zeros(n_int - length(ints), 1)];
                switch node.type
                    case yop.ocp_expr.tp
                        tp = obj.parameterize_timepoint(node, tmp_tp, ...
                            tmp_int, ders);
                        tps = [tps; tp(:)];
                        
                    case yop.ocp_expr.int
                        int = obj.parameterize_integral(node, tmp_tp, ...
                            tmp_int, ders);
                        ints = [ints; int(:)];
                        
                    case yop.ocp_expr.der
                        error(yop.msg.not_implemented);
                        warning(['When implementing dont forget to ' ...
                            'multiply derivative with step length']);
                        
                    otherwise
                        error(yop.msg.unexpected_error);
                end
            end
        end
        
        function val = parameterize_timepoint(obj, tp, tps, ints, ders)
            val = tp.fn(obj.t0, obj.tf, ...
                obj.t.value(tp.timepoint), ...
                obj.x.value(tp.timepoint), ...
                obj.z.value(tp.timepoint), ...
                obj.u.value(tp.timepoint), ...
                obj.p, tps, ints, ders);
        end
        
        function I = parameterize_integral(obj, int, tps, ints, ders)
            I = 0;
            for n=1:obj.N
                yval = [];
                for r=1:length(obj.tau)
                    val_r = int.fn(obj.t0, obj.tf, ...
                        obj.t(n).y(r), ...
                        obj.x(n).y(:, r), ...
                        obj.z(n).evaluate(obj.tau(r)), ...
                        obj.u(n).y, ...
                        obj.p, tps, ints, ders);
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(obj.tau, yval).integrate();
                I = I + lp.evaluate(1)*obj.dt;
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
                if isa(args{k}, 'yop.node')
                    args{k} = obj.value(args{k}, round(mag));
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