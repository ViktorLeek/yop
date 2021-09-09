classdef dc < handle
    properties
        ivals
        points
        polydeg
        
        ocp
        timepoints
        integrals
        
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
        function obj = dc(ocp, ivals, points, polydeg)
            obj.ocp = ocp;
            obj.ivals = ivals;
            obj.points = points;
            obj.polydeg = polydeg;
            obj.timepoints = copy(ocp.timepoints);
            obj.integrals = copy(ocp.integrals);
        end
        
        function obj = build(obj)
            obj.build_vars();
            obj.set_box_bnd();
            obj.discretize_ode();
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
               'x',obj.nx, obj.N+1, obj.d, obj.points);
           
           % Control
           obj.u = cell(obj.N);
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
            for n=1:obj.N
                dx = obj.x(n).differentiate();
                for r=2:obj.d+1
                    xr = obj.x(n).evaluate(obj.tau(r)); % t{n,r}?
                    dxr = dx.evaluate(obj.tau(r));
                    f = obj.ocp.ode(obj.t{n,r}, xr, obj.u{n}, obj.p);
                    obj.add_diffcon(dxr - f);
                end
            end
        end
        
        function obj = add_diffcon(obj, expr)
            if isempty(obj.diffcon)
                obj.diffcon = expr;
            else
                obj.diffcon = [obj.diffcon(:); expr(:)];
            end
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
