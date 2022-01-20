classdef ocp_sol < handle
    properties
        mx_vars
        ids
        t0
        tf
        t
        x
        z
        u
        p
        n_seg
        sol_t
        sol_x
        sol_z
        sol_u
        nx
        nz
        nu
    end
    methods
        function obj = ocp_sol(mx_vars, ids, sol, N, dx, cpx)
            obj.mx_vars = mx_vars;
            obj.ids = ids;
            obj.t0 = sol.t0;
            obj.tf = sol.tf;
            obj.p = sol.p;
            obj.n_seg = sum(N);
            obj.parameterize_polynomials( ...
                sol.t, sol.x, sol.z, sol.u, N, dx, cpx);
            obj.precompute_solution(sol);
            obj.nx = size(sol.x,1);
            obj.nz = size(sol.z,1);
            obj.nu = size(sol.u,1);
        end
        
        
        function obj = precompute_solution(obj, sol)
            % Precompute solution on standard grid
            [zus, uus] = obj.upsample_zu();
            obj.sol_t = sol.t;
            obj.sol_x = sol.x;
            obj.sol_z = zus;
            obj.sol_u = uus;
        end
        
        function [zz,uu] = upsample_zu(obj)
            % Exact upsampling based on solution polynomial
            if isempty(obj.z(1).y)
                zz = [];
            else
                zz = obj.upsample(obj.z);
            end
            if isempty(obj.u(1).y)
                uu = [];
            else
                uu = obj.upsample(obj.u);
            end
        end
        
        function us = upsample(obj, poly)
            rows = size(poly(1).y, 1);
            cols = size(obj.sol_t, 2);
            us = zeros(rows, cols);
            cnt = 1;
            for n=1:obj.n_seg
                for tau = obj.t(n).x
                    us(:,cnt) = poly(n).eval(tau);
                    cnt = cnt + 1;
                end
            end
            us(:,cnt) = poly(n).eval(1);
        end
        
        function obj = parameterize_polynomials(obj, t, x, z, u, N, dx, cpx)
            ip = @(x, y, t0, tf) yop.interpolating_poly(x, y, t0, tf);
            ipe = @(k) yop.interpolating_poly.empty(0, obj.n_seg+k);
            
            taux = yop.ocp_sol.collocation_points(dx, cpx);
            tt=ipe(1); xx=ipe(1); zz=ipe(0); uu=ipe(0);
            cnt=1; ix=1; iz=1; iu=1;
            for r = 1:length(N)
                N_ = N(r);
                dx_ = dx(r);
                taux_ = taux{r};
                tauz_ = taux_(2:end);
                tauu_ = 0;
                for n=1:N_
                    colx = ix : ix + dx_;
                    colz = iz : iz + dx_ - 1;
                    colu = iu;
                    tt0 = t(colx(1));
                    ttf = t(colx(end)+1);
                    tt(cnt) = ip(taux_, t(:, colx), tt0, ttf);
                    xx(cnt) = ip(taux_, x(:, colx), tt0, ttf);
                    zz(cnt) = ip(tauz_, z(:, colz), tt0, ttf);
                    uu(cnt) = ip(tauu_, u(:, colu), tt0, ttf);
                    cnt = cnt + 1;
                    ix = ix + dx_ + 1;
                    iz = iz + dx_;
                    iu = iu + 1;
                end
            end
            assert(cnt==obj.n_seg+1);
            tt(cnt) = ip(0, obj.tf   , obj.tf, obj.tf);
            xx(cnt) = ip(0, x(:, end), obj.tf, obj.tf);
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
            mxv = obj.mx_vars;
            t_id = 0;
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent_initial')
                    mxv{1} = vars{k}.m_value;
                    t_id(end+1) = vars{k}.m_id;
                    
                elseif isa(vars{k}, 'yop.ast_independent_final')
                    mxv{2} = vars{k}.m_value;
                    t_id(end+1) = vars{k}.m_id;
                    
                elseif isa(vars{k}, 'yop.ast_independent')
                    mxv{3} = vars{k}.m_value;
                    t_id(end+1) = vars{k}.m_id;
                    
                end
            end
            
            % Return if variables that are not know are in the ast.
            % This is used in getting values for ocp initial guess, but
            % could also be a user error. Consider fixing in the future.
            IDs = [obj.ids, t_id];
            known_vars = false;
            for k=1:length(vars)
                known_vars = known_vars || any(vars{k}.m_id == IDs);
            end
            if ~known_vars
                v = [];
                return
            end 
            
            % Create functions for the special nodes and expression
            args = { ...
                mxv{:}, ...
                tps.mx_vec(), ...
                ints.mx_vec(), ...
                ders.mx_vec() ...
                };
            
            for node = [tps, ints, ders]
                mx_expr = value(node.ast.m_expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            fn = casadi.Function('fn', args, {value(expr)});
            
            % Compute the numerical values of the special nodes
            [tpv, intv, derv] = obj.comp_sn(sn, ...
                n_elem(tps), n_elem(ints), n_elem(ders));
            
            if isa_reducible(expr)
                v = obj.point(fn, tpv, intv, derv);
            elseif isa_ival(expr)
                v = obj.interval(expr, fn, tpv, intv, derv, mag);
            else
                if mag == 1
                    v = obj.path(fn, tpv, intv, derv);
                else
                    v = obj.pathm(fn, tpv, intv, derv, mag);
                end
            end
        end
        
        function v = point(obj, expr, tps, ints, ders)
            v = full(expr(obj.t0, obj.tf, ...
                obj.t(1).eval(0), ...
                obj.x(1).eval(0), ...
                obj.z(1).eval(0), ...
                obj.u(1).eval(0), ...
                obj.p, tps, ints, ders.evaln(1,0)));
        end
        
        function v = interval(obj, expr, fn, tps, ints, ders, mag)
            [I0, If] = obj.ival_bnds(expr);
            
            % Evaluate at the beginning of the interval
            tt = obj.t.valuenh(I0); 
            xx = obj.x.valuenh(I0);
            zz = obj.z.valuenh(I0);
            uu = obj.u.valuenh(I0);
            dd = ders.valuenh(I0); 
            
            % Evaluate all points within the interval
            [n0, r0, nf, rf] = obj.get_ival_idx(I0, If);
            for n = n0 : min(nf, obj.n_seg)
                tau   = obj.t(n).x;
                r_max = length(tau);
                for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, r_max)
                    tt = [tt, obj.t(n).eval(tau(r))];
                    xx = [xx, obj.x(n).eval(tau(r))];
                    zz = [zz, obj.z(n).eval(tau(r))];
                    uu = [uu, obj.u(n).eval(tau(r))];
                    dd = [dd, ders.evaln(n,tau(r))];
                    if n==nf && r==rf
                        t_nf = obj.t(n).t0;
                        tauf = (If-t_nf)/obj.t(n).dt; % Normalize If
                        tau0 = tau(r);
                        dtau = tauf - tau0;
                    else
                        tau0 = tau(r);
                        if r==r_max
                            tauf = 1;
                        else
                            tauf = obj.t(n).x(r+1);
                        end
                        dtau = tauf - tau0;
                    end
                    for k=1:mag-1 % Magnification
                        tau_k = tau0 + k/mag*dtau;
                        tt = [tt, obj.t(n).eval(tau_k)];
                        xx = [xx, obj.x(n).eval(tau_k)];
                        zz = [zz, obj.z(n).eval(tau_k)];
                        uu = [uu, obj.u(n).eval(tau_k)];
                        dd = [dd, ders.evaln(n,tau_k)];
                    end
                end
            end
            
            % Evaluate last point of interval
            tt = [tt, obj.t.valuenh(If)]; 
            xx = [xx, obj.x.valuenh(If)];
            zz = [zz, obj.z.valuenh(If)];
            uu = [uu, obj.u.valuenh(If)];
            dd = [dd, ders.valuenh(If)]; 
            
            v = full(fn(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
        end
        
        function v = path(obj, expr, tps, ints, ders)
            tt = obj.sol_t; 
            xx = obj.sol_x; 
            zz = obj.sol_z;
            uu = obj.sol_u;
            %dd = [ders.mat, ders(obj.n_seg).eval(1)];
            dd = [ders.mat, ders.evaln(obj.n_seg, 1)];
            v = full(expr(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
        end
        
        function v = pathm(obj, expr, tps, ints, ders, mag)
            cols = length(obj.sol_t)*mag - mag + 1;
            if isempty(ders)
                nd = 0;
            else
                nd = size(ders(1).y,1);
            end
            tt = zeros(     1, cols); 
            xx = zeros(obj.nx, cols); 
            zz = zeros(obj.nz, cols); 
            uu = zeros(obj.nu, cols); 
            dd = zeros(nd , cols);
            cnt = 1;
            for n=1:obj.n_seg
                tau = obj.t(n).x;
                r_max = length(tau);
                for r=1:r_max
                    tt(cnt)   = obj.t(n).eval(tau(r));
                    xx(:,cnt) = obj.x(n).eval(tau(r));
                    zz(:,cnt) = obj.z(n).eval(tau(r));
                    uu(:,cnt) = obj.u(n).eval(tau(r));
                    dd(:,cnt) = ders.evaln(n,tau(r));
                    cnt = cnt + 1;
                    tau0 = tau(r);
                    if r==r_max
                        tauf = 1;
                    else
                        tauf = tau(r+1);
                    end
                    dT = tauf - tau0;
                    for k=1:mag-1 % Magnification
                        tau_k = tau0 + k/mag*dT;
                        tt(cnt)   = obj.t(n).eval(tau_k);
                        xx(:,cnt) = obj.x(n).eval(tau_k);
                        zz(:,cnt) = obj.z(n).eval(tau_k);
                        uu(:,cnt) = obj.u(n).eval(tau_k);
                        dd(:,cnt) = ders.evaln(n,tau_k);
                        cnt = cnt + 1;
                    end
                end
            end
            tt(cnt)   = obj.t(n+1).eval(0);
            xx(:,cnt) = obj.x(n+1).eval(0);
            zz(:,cnt) = obj.z(n).eval(1);
            uu(:,cnt) = obj.u(n).eval(1);
            dd(:,cnt) = ders.evaln(n,1);
            v = full(expr(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
        end
        
%         function v = pathm(obj, expr, tps, ints, ders, mag)
%             tt=[]; xx=[]; zz=[]; uu=[]; dd=[];
%             for n=1:obj.n_seg
%                 tau = obj.t(n).x;
%                 r_max = length(tau);
%                 for r=1:r_max
%                     tt = [tt, obj.t(n).eval(tau(r))];
%                     xx = [xx, obj.x(n).eval(tau(r))];
%                     zz = [zz, obj.z(n).eval(tau(r))];
%                     uu = [uu, obj.u(n).eval(tau(r))];
%                     dd = [dd, ders.evaln(n,tau(r))];
%                     tau0 = tau(r);
%                     if r==r_max
%                         tauf = 1;
%                     else
%                         tauf = tau(r+1);
%                     end
%                     dT = tauf - tau0;
%                     for k=1:mag-1 % Magnification
%                         tau_k = tau0 + k/mag*dT;
%                         tt = [tt, obj.t(n).eval(tau_k)];
%                         xx = [xx, obj.x(n).eval(tau_k)];
%                         zz = [zz, obj.z(n).eval(tau_k)];
%                         uu = [uu, obj.u(n).eval(tau_k)];
%                         dd = [dd, ders.evaln(n,tau_k)];
%                     end
%                 end
%             end
%             tt = [tt, obj.t(n+1).eval(0)];
%             xx = [xx, obj.x(n+1).eval(0)];
%             zz = [zz, obj.z(n).eval(1)];
%             uu = [uu, obj.u(n).eval(1)];
%             dd = [dd, ders.evaln(n,1)];
%             v = full(expr(obj.t0,obj.tf,tt,xx,zz,uu,obj.p,tps,ints,dd));
%         end
        
        function [I0, If] = ival_bnds(obj, expr)
            [I0, If] = get_ival(expr);
            I0 = yop.IF(I0==yop.initial_timepoint, obj.t0, I0);
            I0 = yop.IF(I0==yop.final_timepoint  , obj.tf, I0);
            If = yop.IF(If==yop.initial_timepoint, obj.t0, If);
            If = yop.IF(If==yop.final_timepoint  , obj.tf, If);
        end
        
        function [n0, r0, nf, rf] = get_ival_idx(obj, I0, If)
            
            % First point after I0
            if I0 == obj.tf
                n0 = obj.n_seg + 1;
                r0 = 1;
            else
                for n0 = 1:length(obj.t)
                    tt = obj.t(n0);
                    if tt.t0 <= I0 && I0 < tt.tf
                        tau = (I0 - tt.t0)/tt.dt;
                        taus = tt.x;
                        r0 = find(taus > tau, 1);
                        if isempty(r0)
                            n0 = n0+1;
                            r0 = 1;
                        end
                        break
                    end
                end
            end
            
            % Last point before If
            if If == obj.t0
                nf = 1;
                rf = 1;
                
            else
                for nf = 1:length(obj.t)
                    tt = obj.t(nf);
                    if tt.t0 < If && If <= tt.tf
                        tau = (If - tt.t0)/tt.dt; % tau > 0
                        taus = tt.x;
                        rf = find(taus < tau , 1, 'last');
                        break
                    end
                end
            end
            
        end
        
        function [tps, ints, ders] = comp_sn(obj, sn, n_tp, n_int, n_der)
            tps  = [];
            ints = [];
            ders = obj.init_derivatives(n_der);
            pad = @(dd) [dd;  zeros(n_der  - length(dd), 1)];
            for node = sn
                tmp_tp  = [tps;  zeros(n_tp  - length(tps), 1)];
                tmp_int = [ints; zeros(n_int - length(ints), 1)];
                switch node.type
                    case yop.ocp_expr.tp
                        tp = obj.compute_timepoint( ...
                            node, tmp_tp, tmp_int, ders, pad);
                        tps = [tps; tp(:)];
                        
                    case yop.ocp_expr.int
                        int = obj.compute_integral( ...
                            node, tmp_tp, tmp_int, ders, pad);
                        ints = [ints; int(:)];
                        
                    case yop.ocp_expr.der
                        ders = obj.compute_derivative(node, tmp_tp, ...
                            tmp_int, ders, pad);
                        
                    otherwise
                        error(yop.error.unexpected_error());
                end
            end
        end
        
        function ders = init_derivatives(obj, n_der)
            ders = yop.interpolating_poly.empty(0, obj.n_seg);
            if n_der == 0
                return
            end
            for n=1:obj.n_seg
                tt0 = obj.t(n).t0;
                ttf = obj.t(n).tf;
                tau = obj.t(n).x;
                ders(n) = yop.interpolating_poly(tau, [], tt0, ttf);
            end
        end
        
        function val = compute_timepoint(obj, tp, tps, ints, ders, pad)
            dd = pad(ders.valuenh(tp.timepoint));
            val = full(...
                tp.fn(obj.t0, obj.tf, ...
                obj.t.valuenh(tp.timepoint), ...
                obj.x.valuenh(tp.timepoint), ...
                obj.z.valuenh(tp.timepoint), ...
                obj.u.valuenh(tp.timepoint), ...
                obj.p, tps, ints, dd));
        end
        
        function I = compute_integral(obj, int, tps, ints, ders, pad)
            I = 0;
            for n=1:obj.n_seg
                yval = [];
                tau = obj.t(n).x;
                for r=tau
                    dd = pad(ders.evaln(n,r));
                    val_r = full( ...
                        int.fn(...
                        obj.t0, obj.tf, ...
                        obj.t(n).eval(r), ...
                        obj.x(n).eval(r), ...
                        obj.z(n).eval(r), ...
                        obj.u(n).eval(r), ...
                        obj.p, tps, ints, dd));
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(tau, yval).integrate();
                I = I + lp.evaluate(1)*obj.t(n).dt;
            end
        end
        
        function ders = compute_derivative(obj, der, tps, ints, ders, pad)
            for n=1:obj.n_seg
                yval = [];
                tau = obj.t(n).x;
                for r=tau
                    tt = obj.t(n).eval(r);
                    xx = obj.x(n).eval(r);
                    zz = obj.z(n).eval(r);
                    uu = obj.u(n).eval(r);
                    pp = obj.p;
                    dd = pad(ders(n).eval(r));
                    val_r = full(der.fn( ...
                        obj.t0, obj.tf, tt, xx, zz, uu, pp, tps, ints, dd));
                    yval = [yval, val_r(:)];
                end
                yn = yop.lagrange_polynomial(tau, yval ...
                    ).differentiate().evaluate(tau)/obj.t(n).dt;
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
        function taux = collocation_points(dx, cpx)
            taux = cell(1,length(dx));
            for k=1:length(dx)
                taux{k} = ...
                    full([0, casadi.collocation_points(dx(k), cpx{k})]);
            end
        end
    end

end