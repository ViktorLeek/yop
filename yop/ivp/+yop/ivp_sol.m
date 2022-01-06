classdef ivp_sol < handle
    
    properties
        ivp_vars
        mx_args        
        t
        x
        z
        p
    end
    
    methods
        
        function obj = ivp_sol(ivp_vars, mx_args, t, x, z, p)
            obj.ivp_vars = ivp_vars;
            obj.mx_args = mx_args;
            
            % Numerical solution
            obj.t = t;
            obj.x = x;
            obj.z = z;
            obj.p = p;
        end
        
        
        function v = value(obj, expr)
            
            t_id = 0;
            vars = yop.ivp_sol.find_variables(expr);
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.mx_args{3};
                    t_id = vars{k}.id;
                end
            end
            
            % Return if the variables thats not know are in the ast
            known_vars = false;
            IDs = [obj.ivp_vars.ids, t_id];
            for k=1:length(vars)
                known_vars = known_vars || any(vars{k}.id == IDs);
            end
            if ~known_vars
                v = [];
                return
            end 
            
            obj.ivp_vars.set_mx();
            fn = casadi.Function('fn', obj.mx_args, {fw_eval(expr)});
            
            if isa_reducible(expr)
                v = obj.invariant_value(fn);
            else
                v = obj.variant_value(fn);
            end
            
        end
        
        function v = invariant_value(obj, expr)
            v = full(expr(obj.t(1), obj.t(end), obj.t(1), ...
                obj.x(:,1), obj.z(:,1), obj.p));
        end
        
        function v = variant_value(obj, expr)
            v = full(expr(obj.t(1),obj.t(end),obj.t,obj.x,obj.z,obj.p));
        end
        
        function args = filter(obj, args)
            for k=1:length(args)
                if isa(args{k}, 'yop.ast_node')
                    args{k} = obj.value(args{k});
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
    
    methods (Static)
        function vars = find_variables(expression)
            % Find and add all special nodes of the expression
            [tsort, N] = topological_sort(expression);
            
            % Find special nodes, maintain topological order
            vars = {};
            for n=1:N                
                if isa(tsort{n}, 'yop.ast_variable')
                    vars{end+1} = tsort{n};
                end
            end
        end
    end
    
end