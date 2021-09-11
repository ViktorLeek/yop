classdef dc < handle
    properties
        ivals
        points
        polydeg
        
        ocp
        timepoints
        integrals
        eq
        ieq
        objective
        
        tau
        t0
        tf
        t
        x
        u
        p
        tp
        int
        
        t0_ub
        t0_lb
        tf_ub
        tf_lb
        x_ub
        x_lb
        u_ub
        u_lb
        p_ub
        p_lb
        
        diffcon
        
    end
    methods
        function obj = dc(ocp, ivals, polydeg, points)
            obj.ocp = ocp;
            obj.ivals = ivals;
            obj.points = points;
            obj.polydeg = polydeg;
            obj.timepoints = copy(ocp.timepoints);
            obj.integrals = copy(ocp.integrals);
            obj.eq = copy(ocp.eq);
            obj.ieq = copy(ocp.ieq);
            obj.objective = copy(ocp.objective);
        end
        
        function nlp = build(obj)
            obj.build_vars();
            obj.set_box_bnd();
            obj.discretize_ode();
            obj.parameterize_timepoints_and_integrals();
            obj.parameterize_expressions();
            d = obj.diffcon;
            g = obj.eq.vertcat_disc();
            h = obj.ieq.vertcat_disc();
            nlp.f = obj.objective.disc;
            nlp.x = obj.w;
            nlp.g = [d; g; h];
            nlp.x_lb = obj.w_lb;
            nlp.x_ub = obj.w_ub;
            nlp.g_ub = [zeros(size([d;g])); zeros(size(h))];
            nlp.g_lb = [zeros(size([d;g])); -inf(size(h))];
        end
        
        function obj = build_vars(obj)
           obj.tau = [0, ...
               casadi.collocation_points(obj.polydeg, obj.points)];
           
           % Independent
           obj.t0 = casadi.MX.sym('t0');
           obj.tf = casadi.MX.sym('tf');
           obj.t = cell(obj.N+1, obj.d+1);
           for n=1:obj.N
               for r=1:(obj.d+1)
                   obj.t{n,r} = obj.t0 + obj.h*(n-1) + obj.h*obj.tau(r);
               end
           end
           obj.t{obj.N+1,1} = obj.tf;
           
           % State
           obj.x = yop.collocated_state( ...
               'x', obj.nx, obj.N, obj.d, obj.points);
           
           % Control
           obj.u = cell(obj.N, 1);
           for n=1:obj.N
               obj.u{n} = casadi.MX.sym(['u_', num2str(n)], obj.nu);
           end
           
           % Parameter
           obj.p = casadi.MX.sym('p', obj.np);
           
           % Timepoint
           obj.tp = casadi.MX.sym('tp', obj.ntp);
           
           % Integral
           obj.int = casadi.MX.sym('int', obj.nint);
           
        end
        
        function obj = set_box_bnd(obj)
            obj.t0_lb = obj.ocp.t0_lb;
            obj.t0_ub = obj.ocp.t0_ub;
            obj.tf_lb = obj.ocp.tf_lb;
            obj.tf_ub = obj.ocp.tf_ub;
            
            reps = obj.N*(obj.d + 1) + 1;
            obj.x_lb = repmat(obj.ocp.x_lb, reps, 1);
            obj.x_ub = repmat(obj.ocp.x_ub, reps, 1);
            obj.x_lb(1:obj.nx) = obj.ocp.x0_lb;
            obj.x_ub(1:obj.nx) = obj.ocp.x0_ub;
            obj.x_lb(end-obj.nx+1:end) = obj.ocp.xf_lb;
            obj.x_ub(end-obj.nx+1:end) = obj.ocp.xf_ub;
            
            obj.u_lb = repmat(obj.ocp.u_lb, obj.N, 1);
            obj.u_ub = repmat(obj.ocp.u_ub, obj.N, 1);
            obj.u_lb(1:obj.nu) = obj.ocp.u0_lb;
            obj.u_ub(1:obj.nu) = obj.ocp.u0_ub;
            obj.u_lb(end-obj.nu+1:end) = obj.ocp.uf_lb;
            obj.u_ub(end-obj.nu+1:end) = obj.ocp.uf_ub;
                        
            obj.p_lb = obj.ocp.p_lb;
            obj.p_ub = obj.ocp.p_ub;
        end
        
        function obj = discretize_ode(obj)
            %             for n=1:obj.N
            %                 xx = obj.x(n).evaluate(obj.tau);
            %                 dx = obj.x(n).differentiate().evaluate(obj.tau);
            %                 for r=2:obj.d+1
            %                     f = obj.ocp.ode(obj.t{n,r}, xx(r), obj.u{n}, obj.p);
            %                     obj.add_diffcon(dx(r) - obj.h*f);
            %                 end
            %             end
            for n=1:obj.N
                dx = obj.x(n).differentiate();
                for r=2:obj.d+1
                    xr = obj.x(n).evaluate(obj.tau(r)); % t{n,r}?
                    dxr = dx.evaluate(obj.tau(r));
                    f = obj.ocp.ode(obj.t{n,r}, xr, obj.u{n}, obj.p);
                    obj.add_diffcon(dxr - obj.h*f);
                end
            end
            
            for n=1:obj.N
                obj.add_diffcon( ...
                    obj.x(n).evaluate(1) - obj.x(n+1).evaluate(0));
            end
            
        end
        
        function obj = parameterize_timepoints_and_integrals(obj)
            for tp_k = obj.timepoints
                [tt,xx,uu,pp] = obj.vars_at(tp_k.timepoint);
                tvec = obj.timepoints.get_disc();
                ivec = obj.integrals.get_disc();
                tmp = tp_k.fn(tt,xx,uu,pp,tvec, ivec);
                tp_k.disc = tmp(:);
            end
            tvec = obj.timepoints.get_disc();
            test = casadi.Function('tvec', {obj.w}, {tvec});
            
            for int_k = obj.integrals
                ivec = obj.integrals.get_disc();
                I = 0;
                for n=1:obj.N
                    val = [];
                    uu = obj.u{n};
                    pp = obj.p;
                    for r=1:obj.d+1
                        tt = obj.t{n,r};
                        xx = obj.x(n).evaluate(obj.tau(r));
                        val_r = int_k.fn(tt,xx,uu,pp,tvec,ivec);
                        val = [val, val_r(:)];
                    end
                    lp = yop.lagrange_polynomial(obj.tau, val).integrate();
                    I = I + lp.evaluate(1);
                end
                int_k.disc = I(:);
            end
        end
        
        function obj = parameterize_expressions(obj)
            tvec = obj.timepoints.get_disc();
            ivec = obj.integrals.get_disc();
            for c = [obj.objective, obj.eq, obj.ieq]
                if all(is_transcription_invariant(c))
                    tt = obj.t{1,1};
                    xx = obj.x(1).evaluate(0);
                    uu = obj.u{1};
                    pp = obj.p;
                    tmp = c.fn(tt,xx,uu,pp,tvec,ivec);
                    c.disc = tmp(:);
                else
                    % IMPLEMENT HARD
                    grid_val = cell(obj.N+1, 1);
                    for n=1:obj.N
                        tt = obj.t{n,1};
                        xx = obj.x(n).evaluate(0);
                        uu = obj.u{n};
                        pp = obj.p;
                        tmp = c.fn(tt,xx,uu,pp,tvec,ivec);
                        grid_val{n} = tmp(:);
                    end
                    tt = obj.t{n+1,1};
                    xx = obj.x(n+1).evaluate(0);
                    uu = obj.u{n};
                    pp = obj.p;
                    tmp = c.fn(tt,xx,uu,pp,tvec,ivec);
                    grid_val{n+1} = tmp(:);
                    c.disc = vertcat(grid_val{:});
                end
            end
        end
        
        function [t, x, u, p] = vars_at(obj, tp)
                
            if isa(tp, 'yop.ast_independent_initial')
                t = obj.t0;
                x = obj.x(1).y(:,1); % Need an abstraction here
                u = obj.u{1};
                p = obj.p;
                return;
            elseif isa(tp, 'yop.ast_independent_final')
                t = obj.tf;
                x = obj.x(end).y; % y-value at N+1 is a vector 
                u = obj.u{end};
                p = obj.p;
                return;
            end
            
            [fixed, T] = fixed_horizon(obj.ocp);            
            
            if ~fixed
                error('[Yop] Problem does not have a fixed horizon')
            end
            
            dt = T/obj.N; % grid step_length    
            if yop.EQ(rem(tp,dt), 0, 1e-3) % 1e-3 should be a setting (also for DMS)
                % grid point
                n = round(tp/dt);
                t = obj.t{n};
                x = obj.x(n).y(:,1);
                u = obj.u{min(obj.N, n)};
                p = obj.p;
                
            else
                % With the first statement and this floor operation n!=N+1,
                % unless t is outside the horizon, in which case an error
                % is desirable.
                % Integrate up to timepoint
                n = floor(tp/dt);
                x = obj.x(n).evaluate(tp);
                u = obj.u{n};
                p = obj.p;
            end
        end
        
        function obj = add_diffcon(obj, expr)
            if isempty(obj.diffcon)
                obj.diffcon = expr;
            else
                obj.diffcon = [obj.diffcon(:); expr(:)];
            end
        end
        
        function v = w(obj)
            xx = [];
            for xk=obj.x
                xx = [xx(:); xk.y(:)];
            end
            v = vertcat(obj.t0, obj.tf, xx, obj.u{:}, obj.p);
        end
        
        function bd = w_lb(obj)
            bd = vertcat(obj.t0_lb, obj.tf_lb, obj.x_lb, obj.u_lb, obj.p_lb);
        end
        
        function bd = w_ub(obj)
            bd = vertcat(obj.t0_ub, obj.tf_ub, obj.x_ub, obj.u_ub, obj.p_ub);
        end
        
        function n = N(obj)
            n = obj.ivals;
        end
        
        function n = d(obj)
            n = obj.polydeg;
        end
        
        function n = h(obj)
            n = (obj.tf-obj.t0)/obj.N;
        end
        
        function n = nx(obj)
            n = obj.ocp.n_x;
        end
        
        function n = nu(obj)
            n = obj.ocp.n_u;
        end
        
        function n = np(obj)
            n = obj.ocp.n_p;
        end
        
        function n = ntp(obj)
            n = obj.ocp.n_tp;
        end
        
        function n = nint(obj)
            n = obj.ocp.n_int;
        end
    end
end
