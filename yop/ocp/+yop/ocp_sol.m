classdef ocp_sol < handle
    properties
        sol
        solver
        nlp
        
        ocp_t
        ocp_x
        ocp_z
        ocp_u
        ocp_p
    end
    methods
        function obj = ocp_sol(sol, solver, nlp, t, x, z, u, p)
            obj.sol = sol;
            obj.solver = solver;
            obj.nlp = nlp;
            obj.ocp_t = t;
            obj.ocp_x = x;
            obj.ocp_z = z;
            obj.ocp_u = u;
            obj.ocp_p = p;
        end
        
        function obj = save(obj, filename)
            res = struct;
            
            % This part is sufficient to reconstruct the results
            res.w = full(obj.sol.x);
            res.N = obj.nlp.N;
            res.d = obj.nlp.d;
            res.cp = obj.nlp.cp;
            
            % But as ordering of the states matter for interpreting the
            % results, the following convenience entries are added, they
            % also provide means of semantically meaningful symbolic
            % plotting.
            x_name = {};
            for x=obj.ocp_x
                x_name{end+1} = x.var.name;
            end
            
            z_name = {};
            for z=obj.ocp_z
                z_name{end+1} = z.var.name;
            end
            
            u_name = {};
            for u=obj.ocp_u
                u_name{end+1} = u.var.name;
            end
            
            p_name = {};
            for p=obj.ocp_p
                p_name{end+1} = p.var.name;
            end
            
            res.x_name = x_name;
            res.z_name = z_name;
            res.u_name = u_name;
            res.p_name = p_name;
            
            time = casadi.Function('t', {obj.nlp.w}, {mat(obj.nlp.t)});
            state = casadi.Function('x', {obj.nlp.w}, {mat(obj.nlp.x)});
            algebraic = casadi.Function('x', {obj.nlp.w}, {[]});
            control = casadi.Function('u', {obj.nlp.w}, {mat(obj.nlp.u)});
            parameter = casadi.Function('p', {obj.nlp.w}, {obj.nlp.p});
            
            res.t_sol = full(time(obj.sol.x));
            res.x_sol = full(state(obj.sol.x));
            res.z_sol = full(algebraic(obj.sol.x));
            res.u_sol = full(control(obj.sol.x));
            res.p_sol = full(parameter(obj.sol.x));
            
            res.nx = n_elem(obj.ocp_x);
            res.nu = n_elem(obj.ocp_u);
            res.nz = 0;
            res.np = n_elem(obj.ocp_p);
            
            if strcmp(filename(end-3:end), '.mat')
                save(filename, '-struct', 'res');
            else
                save([filename, '.mat'], '-struct', 'res');
            end
        end
        
        function obj = reinit(obj,t,x,z,u,p)
        end
        
        function v = value(obj, expr)
            % This discretizes the expression in the same way that the
            % transcription method does it. The idea is to see what the
            % solver sees.
            
            expr = yop.ocp_expr(expr);
            
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr.ast);
            
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.ocp_t.mx;
                end
            end
            
            args = { ...
                mx_vec(obj.ocp_t), ...
                mx_vec(obj.ocp_x), ...
                mx_vec(obj.ocp_u), ...
                mx_vec(obj.ocp_p), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx([obj.ocp_t, obj.ocp_x, obj.ocp_u, obj.ocp_p]);
            set_mx([tps, ints, ders]);
            
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            
            expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});
            
            t0 = full(obj.sol.x(1));
            tf = full(obj.sol.x(2));
            [tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), ...
                n_elem(ints), obj.nlp.N, obj.nlp.tau, obj.nlp.dt, t0, ...
                tf, obj.nlp.t, obj.nlp.x, obj.nlp.u, obj.nlp.p);
            
            disc = yop.parameterize_expression(expr, obj.nlp.N, ...
                obj.nlp.tau, obj.nlp.t, obj.nlp.x, obj.nlp.u, ...
                obj.nlp.p, tpv, intv, []);
            
            f = casadi.Function('f', {obj.nlp.w}, {disc});
            v = full(f(obj.sol.x));
        end
        
        function v = valuex(obj, expr, refinement)
            expr = yop.ocp_expr(expr);
            
            %% This could be part of the contructor and tt, xx ... could be
            %  instance variables
            time = casadi.Function('t', {obj.nlp.w}, {mat(obj.nlp.t)});
            tt = full(time(obj.sol.x));
            
            state = casadi.Function('x', {obj.nlp.w}, {mat(obj.nlp.x)});
            xx = full(state(obj.sol.x));
            
            control = casadi.Function('u', {obj.nlp.w}, {mat(obj.nlp.u)});
            uu = full(control(obj.sol.x));
            
            parameter = casadi.Function('p', {obj.nlp.w}, {obj.nlp.p});
            pp = full(parameter(obj.sol.x));
            
            tlp = yop.collocated_time(tt(1), tt(end), obj.nlp.N);
            
            y = cell(obj.nlp.N+1,1);
            for n=1:obj.nlp.N
                y{n} = xx(:, (obj.nlp.d+1)*n-obj.nlp.d:(obj.nlp.d+1)*n);
            end
            y{obj.nlp.N+1} = xx(:,end);
            xlp = yop.collocated_expression(obj.nlp.N, obj.nlp.tau, y);
            
            y = cell(obj.nlp.N+1,1);
            for n=1:obj.nlp.N
                y{n} = uu(:, n);
            end
            y{obj.nlp.N+1} = uu(:,end);
            ulp = yop.collocated_expression(obj.nlp.N, 0, y);
            %%
            
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr.ast);
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.ocp_t.mx;
                end
            end
            
            args = { ...
                mx_vec(obj.ocp_t), ...
                mx_vec(obj.ocp_x), ...
                mx_vec(obj.ocp_u), ...
                mx_vec(obj.ocp_p), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx([obj.ocp_t, obj.ocp_x, obj.ocp_u, obj.ocp_p]);
            set_mx([tps, ints, ders]);
            
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            
            expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});
            
            t0 = full(obj.sol.x(1));
            tf = full(obj.sol.x(2));
            [tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), ...
                n_elem(ints), obj.nlp.N, obj.nlp.tau, obj.nlp.dt, t0, ...
                tf, obj.nlp.t, obj.nlp.x, obj.nlp.u, obj.nlp.p);
            
            tpf = casadi.Function('f', {obj.nlp.w}, {tpv});
            intf = casadi.Function('f', {obj.nlp.w}, {intv});
            
            TP = tpf(obj.sol.x);
            I = intf(obj.sol.x);
            
            if expr.is_transcription_invariant
                v = expr.fn( ...
                    tlp(1).evaluate(0), ...
                    xlp(1).evaluate(0), ...
                    ulp(1).evaluate(0), ...
                    pp, TP, I, []);
            else
                t = tt(1:2:end);
                dt = t(2)-t(1);
                tvec = t;
                for k=1:refinement-1
                    tvec = [tvec; t+k*dt/refinement];
                end
                v = [];
                for tp=tvec(:)'
                    vk = expr.fn( ...
                        tlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        xlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        ulp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        pp, TP, I, []);
                    v = [v, vk];
                end
            end
            v = full(v);
        end
        
        function args = filter(obj, input)
            args = {};
            refine = false;
            refinement = 0;
            k = 1;
            while k <= length(input)
                if strcmp(input{k}, 'refine')
                    refine = true;
                    refinement = input{k+1};
                    k = k + 2;
                else
                    args{end+1} = input{k};
                    k = k + 1;
                end
            end
            
            for k=1:length(args)
                if isa(args{k}, 'yop.node')
                    if refine
                        args{k} = obj.valuex(input{k}, refinement);
                    else
                        args{k} = obj.value(input{k});
                    end
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