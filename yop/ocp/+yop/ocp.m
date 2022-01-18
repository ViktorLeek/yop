classdef ocp < handle
    properties
        name
        
        % Variables
        independent
        independent0
        independentf
        states
        algebraics
        controls
        parameters
        
        % Objective
        objective
        
        % Dynamics
        ode
        alg
        
        % Constraints
        path
        path_hard
        path_ival
        point
        
        % User input
        ode_eqs
        alg_eqs
        ec_eqs
        ec_hard_eqs
        ec_ival_eqs
        ec_point_eqs
        iec_eqs
        iec_hard_eqs
        iec_ival_eqs
        iec_point_eqs
        guess
        
    end
    
    properties (Hidden) % For parsing
        snodes % Topological sort of special nodes
        tps  % Timepoints
        ints % Integrals
        ders % Derivativtes (time varying value)
        visited = [] % special nodes that has been visited - speedup
    end
    
    properties %(Hidden) % Multiphase
        m_id     
        m_value  
        m_phases 
        m_parent 
        m_nlp
    end
    
    %% User interface
    methods 
        function obj = ocp(name)
            if nargin == 1
                obj.name = name;
            else
                obj.name = 'Optimal Control Problem';
            end
            obj.ode_eqs = {};
            obj.alg_eqs = {};
            obj.ec_eqs  = {};
            obj.ec_hard_eqs  = {};
            obj.ec_ival_eqs  = {};
            obj.ec_point_eqs = {};
            obj.iec_eqs = {};
            obj.iec_hard_eqs  = {};
            obj.iec_ival_eqs  = {};
            obj.iec_point_eqs = {};
            
            obj.snodes = yop.ocp_expr.empty(1,0);
            obj.tps    = yop.ocp_expr.empty(1,0);
            obj.ints   = yop.ocp_expr.empty(1,0);
            obj.ders   = yop.ocp_expr.empty(1,0);
            
            obj.independent  = yop.ocp_var.empty(1,0);
            obj.independent0 = yop.ocp_var.empty(1,0);
            obj.independentf = yop.ocp_var.empty(1,0);
            obj.states       = yop.ocp_var.empty(1,0);
            obj.algebraics   = yop.ocp_var.empty(1,0);
            obj.controls     = yop.ocp_var.empty(1,0);
            obj.parameters   = yop.ocp_var.empty(1,0);
            
            obj.m_id = yop.ocp.get_id();
            obj.m_value = yop.cx(['ocp', num2str(obj.m_id)]);
            obj.m_phases = obj;
            obj.m_parent = obj; % Self means no parent
            
        end
        
        function obj = min(obj, expr)
            obj.set_objective(expr);
        end
        
        function obj = max(obj, expr)
            obj.set_objective(-expr);
        end
        
        function obj = st(obj, varargin)
            for k=1:length(varargin)
                obj.parse_constraint(varargin{k});
            end
        end
        
        function obj = hard(obj, varargin)
            for k=1:length(varargin)
                varargin{k} = hard(varargin{k});
            end
            obj.st(varargin{:});
        end
        
        function sol = solve(obj, varargin)
            ip = inputParser();
            ip.FunctionName = "yop.ocp/solve";
            ip.addParameter('intervals', yop.defaults.control_invervals);
            ip.addParameter('degree', yop.defaults.polynomial_degree);
            ip.addParameter('points', yop.defaults.collocation_points);
            ip.addParameter('guess', []);
            ip.addParameter('solver', 'ipopt');
            ip.addParameter('opts', struct());
            ip.parse(varargin{:});
            
            N  = ip.Results.intervals;
            d  = ip.Results.degree;
            cp = ip.Results.points;
            opts = ip.Results.opts;
            solver = ip.Results.solver;
            obj.guess = ip.Results.guess;
            
            obj.to_canonical();
            nlp = yop.direct_collocation(obj, N, d, cp);
            
            solver = casadi.nlpsol('solver', solver, ...
                struct('f', nlp.J, 'x', nlp.w, 'g', nlp.g), opts);
            nlp_sol = solver( ...
                'x0', nlp.w0, ...
                'lbx', nlp.w_lb, ...
                'ubx', nlp.w_ub, ...
                'ubg', nlp.g_ub, ...
                'lbg', nlp.g_lb ...
                );
            
            w_opt = struct;
            w_opt.t0 = obj.descale_t0(nlp.ocp_t0(nlp_sol.x));
            w_opt.tf = obj.descale_tf(nlp.ocp_tf(nlp_sol.x));
            w_opt.t  = obj.descale_t (nlp.ocp_t (nlp_sol.x));
            w_opt.x  = obj.descale_x (nlp.ocp_x (nlp_sol.x));
            w_opt.z  = obj.descale_z (nlp.ocp_z (nlp_sol.x));
            w_opt.u  = obj.descale_u (nlp.ocp_u (nlp_sol.x));
            w_opt.p  = obj.descale_p (nlp.ocp_p (nlp_sol.x));
            
            sol = yop.ocp_sol(obj.mx_vars(), obj.ids, w_opt, N, d, {cp});
            yop.progress.ocp_solved(solver.stats.success);
        end
    end
    
    %% Multiphase
    methods
        function varargout = solve_multiphase(obj, N, d, cp, solver, opts)
            
            % Build individual NLPs
            
            obj.build_nlps(N, d, cp);
            nlp = obj.merge_nlps();
            nlp_sol = solve_nlp(solver, opts);            
            obj.nlp_to_ocp_maps(nlp);
            %osol = obj.ocp_sol(nlp_sol);
            psol = obj.phase_sol(nlp_sol);
            
            varargout = {osol};
            %             for
            % Merge variables
        end
        
        function sol = ocp_sol(obj, nlp_sol)
            % Solution to the complete OCP
            t=[]; x=[]; z=[]; u=[]; p=[];
            for p = obj.m_phases
                tt = p.descale_t( p.m_nlp.ocp_t(nlp_sol.x) );
                xx = p.descale_x( p.m_nlp.ocp_x(nlp_sol.x) );
                zz = p.descale_z( p.m_nlp.ocp_z(nlp_sol.x) );
                uu = p.descale_u( p.m_nlp.ocp_u(nlp_sol.x) );
                t = [t, tt(2:end)];
                x = [x, xx(:, 2:end)];
                z = [z, zz];
                u = [u, uu];
            end
            p1 = obj.m_phases(1);
            pf = obj.m_phases(end);
            t0 = p1.descale_t0(p1.m_nlp.ocp_t0(nlp_sol.x));
            tf = pf.descale_tf(pf.m_nlp.ocp_tf(nlp_sol.x));
            sol = struct('t0',t0,'tf',tf,'t',t,'x',x,'z',z,'u',u,'p',p);
        end
        
        function s = phase_sol(obj, nlp_sol)
            s = struct;
            s.t0=[]; s.tf=[]; s.t=[]; s.x=[]; s.z=[]; s.u=[]; s.p=[];
            for p=obj.m_phases
                s(end+1).t0 = p.descale_t0(p.m_nlp.ocp_t0(nlp_sol.x));
                s(end+1).tf = p.descale_tf(p.m_nlp.ocp_tf(nlp_sol.x));
                s(end+1).t  = p.descale_t (p.m_nlp.ocp_t (nlp_sol.x));
                s(end+1).x  = p.descale_x (p.m_nlp.ocp_x (nlp_sol.x));
                s(end+1).z  = p.descale_z (p.m_nlp.ocp_z (nlp_sol.x));
                s(end+1).u  = p.descale_u (p.m_nlp.ocp_u (nlp_sol.x));
                s(end+1).p  = p.descale_p (p.m_nlp.ocp_p (nlp_sol.x));
            end
        end
        
        function build_nlps(obj, N, d, cp)
            for k=1:length(obj.m_phases)
                obj.m_phases(k).to_canonical();
                obj.m_phases(k).m_nlp = ...
                    yop.direct_collocation(N(k), d(k), cp{k});
            end
        end
        
        function nlp_sol = solve_nlp(~, solver, opts)
            prob = struct('f', nlp.J, 'x', nlp.w, 'g', nlp.g);
            solver = casadi.nlpsol('solver', solver, prob, opts);
            nlp_sol = solver( ...
                'x0' , nlp.w0, ...
                'lbx', nlp.w_lb, ...
                'ubx', nlp.w_ub, ...
                'ubg', nlp.g_ub, ...
                'lbg', nlp.g_lb ...
                );
        end
        
        function nlp = merge_nlps(obj)
            parent = obj.parent_position;
            nlp    = yop.nlp;
            nlp.J  = obj.nlp_objective();
            for k=1:length(obj.m_phases)
                pk = obj.m_phases(k);
                nlp.w    = [nlp.w   ; pk.m_nlp.w   ];
                nlp.w0   = [nlp.w0  ; pk.m_nlp.w0  ];
                nlp.w_ub = [nlp.w_ub; pk.m_nlp.w_ub];
                nlp.w_lb = [nlp.w_lb; pk.m_nlp.w_lb];
                nlp.g    = [nlp.g   ; pk.m_nlp.g   ];
                nlp.g_ub = [nlp.g_ub; pk.m_nlp.g_ub];
                nlp.g_lb = [nlp.g_lb; pk.m_nlp.g_lb];
                obj.continuity(nlp, pk, obj.m_phases(parent(k)));
            end
        end
        
        function J = nlp_objective(obj)
            % First create a function object for the multiphase objective
            args = arrayfun( ...
                @(e) e.m_value, J.m_phases, 'UniformOutput', false);
            Jfn = casadi.Function('J', args, {obj.m_value});
            
            % Evaluate the function using the individual nlp objectives
            Js = arrayfun( ...
                @(e) e.m_nlp.J, J.m_phases, 'UniformOutput', false);
            J = Jfn(Js{:});
            
        end
        
        function continuity(~, nlp, parent, child)
            if parent.m_id == child.m_id
                % Circular OCPs formualted as multi-phase problems are not
                % supported
                return;
            end
            
            nlp_p = parent.m_nlp;
            nlp_c = child.m_nlp;
            
            time_pos    = nlp_c.t0 - nlp_c.tf; % constratint <= 0
            time_cont   = nlp_p.tf - nlp_c.t0; % constratint == 0
            state_cont  = nlp_p.x(end).eval(0) - nlp_c.x(1).eval(0);
            param_const = nlp_p.p - nlp_c.p;
            
            eq   = [time_cont; state_cont; param_const];
            ieq  = time_pos;
            g    = [eq; ieq];
            g_ub = [zeros(size(eq)); zeros(size(ieq))];
            g_lb = [zeros(size(eq));  -inf(size(ieq))];
            
            nlp.g    = [nlp.g   ; g   ];
            nlp.g_ub = [nlp.g_ub; g_ub];
            nlp.g_lb = [nlp.g_lb; g_lb];
        end
        
        function parent = parent_position(obj)
            % convert object to its id
            parent = arrayfun(@(e) e.m_id, obj.m_parent);
            for k=1:length(obj.m_phases)
                % IDs are positive so positions are mapped to negative
                % values
                % Replace parents id with its position in the phase list,
                % but with a negative value so as to avoid mixing indices
                % and ids.
                parent(parent==obj.m_phases(k).m_id) = -k;
            end
            % Convert to positive indices
            parent = abs(parent);
        end
        
        function nlp_to_ocp_maps(~, nlp)
            for p=1:obj.m_phases
                % Function that maps from multiphase nlp variable vector to
                % single phase variable vector
                f = casadi.Function('f', {nlp.w}, {p.m_nlp.w});
                
                % Functions that map from nlp vector to ocp variables
                p.m_nlp.ocp_t0 = @(w) p.m_nlp.ocp_t0(full(f(w)));
                p.m_nlp.ocp_tf = @(w) p.m_nlp.ocp_tf(full(f(w)));
                p.m_nlp.ocp_t  = @(w) p.m_nlp.ocp_t (full(f(w)));
                p.m_nlp.ocp_x  = @(w) p.m_nlp.ocp_x (full(f(w)));
                p.m_nlp.ocp_z  = @(w) p.m_nlp.ocp_z (full(f(w)));
                p.m_nlp.ocp_u  = @(w) p.m_nlp.ocp_u (full(f(w)));
                p.m_nlp.ocp_p  = @(w) p.m_nlp.ocp_p (full(f(w)));
            end
        end
        
        
        function obj = set_guess(obj, guess)
           obj.guess = guess; 
        end
    end
    
    %% Parsing
    methods 
        
        function set_objective(obj, expr)
            % Topological sort of expression in order to find variables,
            % timepoints, integrals and derivatives.
            tp_int = obj.find_special_nodes(expr);
            
            % Error handling - Test if timedependent variables reach the
            % objective without being in a timepoint or integral. By
            % making a second topological sort, where timepoints and
            % integrals have been removed, it is tested whether time
            % dependent variabels are reached in the ast. If a time
            % dependent variable is inside a timepoint or integral
            % it will not be visited by the topological sort.
            if yop.settings.errors
                visited_ = yop.get_ids(tp_int);
                [sort, N] = topological_sort(expr, visited_);
                err_nodes = {};
                for n=1:N
                    time_varying = ...
                        isa(sort{n}, 'yop.ast_state') || ...
                        isa(sort{n}, 'yop.ast_control') || ...
                        isa(sort{n}, 'yop.ast_algebraic') || ...
                        isa(sort{n}, 'yop.ast_independent');
                    if time_varying
                        err_nodes{end+1} = sort{n};
                    end
                end
                
                if ~isempty(err_nodes)
                    error(yop.error.timevarying_objective(err_nodes));
                end
            end
            
            % Passed error check - assign value
            obj.objective.ast = expr;
        end
        
        function parse_constraint(obj, c)
            obj.find_special_nodes(c);
            ssr = yop.to_ssr(c);
            for k=1:length(ssr)
                obj.classify_constraint(ssr{k});
            end
        end
        
        function classify_constraint(obj, ssr)
            
            lhs = ssr.m_lhs;
            rhs = ssr.m_rhs;
            t0 = yop.initial_timepoint();
            tf = yop.final_timepoint();
            
            isa_eq = isa(ssr, 'yop.ast_eq');
            isa_le = isa(ssr, 'yop.ast_le') || isa(ssr, 'yop.ast_lt');
            isa_ge = isa(ssr, 'yop.ast_ge') || isa(ssr, 'yop.ast_gt');
            invariant = isa_reducible(ssr);
            hard = is_hard(ssr);
            
            isnum_lhs = isa_numeric(lhs);
            num_lhs   = numval(lhs);
            ival_lhs  = isa_ival(lhs);
            istp_lhs  = isa_timepoint(lhs);
            der_lhs   = isa_der(lhs);
            var_lhs   = isa_variable(lhs);
            fnh_lhs   = isa(ssr.m_lhs, 'function_handle');
            [type_lhs, id_lhs] = Type(lhs);
            state_lhs   = type_lhs == yop.var_type.state;
            control_lhs = type_lhs == yop.var_type.control;
            
            isnum_rhs = isa_numeric(rhs);
            num_rhs   = numval(rhs);
            ival_rhs  = isa_ival(rhs);
            istp_rhs  = isa_timepoint(rhs);
            der_rhs   = isa_der(rhs);
            var_rhs   = isa_variable(rhs);
            fnh_rhs   = isa(ssr.m_rhs, 'function_handle');
            [type_rhs, id_rhs] = Type(rhs);
            state_rhs   = type_rhs == yop.var_type.state;
            control_rhs = type_rhs == yop.var_type.control;
            
            if (ival_lhs || ival_rhs) && isa_eq
                % interval equality constraint
                obj.ec_ival_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif (isa_ival(lhs) || isa_ival(rhs))
                % interval inequality constraint
                obj.iec_ival_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif der_lhs && state_lhs && isa_eq
                % der(x) == expr
                obj.ode_eqs{end+1} = ssr;
                obj.remove_state_der(lhs.m_der);
                
            elseif der_rhs && state_rhs && isa_eq
                % expr == der(x)
                c = get_constructor(ssr);
                obj.ode_eqs{end+1} = c(ssr.m_rhs, ssr.m_lhs);
                obj.remove_state_der(rhs.m_der);
                
            elseif isa_alg(ssr) && isa_eq
                % alg(expr1 == expr2)
                obj.alg_eqs{end+1} = ssr;
                
            elseif istp_lhs && lhs.m_t0==t0 && isnum_rhs && isa_eq && (state_lhs || control_lhs)
                % v(t0) == num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ub0 = bnd;
                var.lb0 = bnd;
                
            elseif istp_lhs && lhs.m_t0==t0 && isnum_rhs && isa_le && (state_lhs || control_lhs)
                % v(t0) <= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ub0 = bnd;
                
            elseif istp_lhs && lhs.m_t0==t0 && isnum_rhs && isa_ge && (state_lhs || control_lhs)
                % v(t0) >= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.lb0 = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_t0==t0 && isa_eq && (state_rhs || control_rhs)
                % num == v(t0)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ub0 = bnd;
                var.lb0 = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_t0==t0 && isa_le && (state_rhs || control_rhs)
                % num <= v(t0)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.lb0 = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_t0==t0 && isa_ge && (state_rhs || control_rhs)
                % num >= v(t0)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ub0 = bnd;
                
            elseif istp_lhs && lhs.m_tf==tf && isnum_rhs && isa_eq && (state_lhs || control_lhs)
                % v(tf) == num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ubf = bnd;
                var.lbf = bnd;
                
            elseif istp_lhs && lhs.m_tf==tf && isnum_rhs && isa_le && (state_lhs || control_lhs)
                % v(tf) <= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ubf = bnd;
                
            elseif istp_lhs && lhs.m_tf==tf && isnum_rhs && isa_ge && (state_lhs || control_lhs)
                % v(tf) >= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.lbf = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_tf==tf && isa_eq && (state_rhs || control_rhs)
                % num == v(tf)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ubf = bnd;
                var.lbf = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_tf==tf && isa_le && (state_rhs || control_rhs)
                % num <= v(tf)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.lbf = bnd;
                
            elseif isnum_lhs && istp_rhs && rhs.m_tf==tf && isa_ge && (state_rhs || control_rhs)
                % num >= v(tf)
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ubf = bnd;
                
            elseif var_lhs && isnum_rhs && isa_eq
                % v == num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ub = bnd;
                var.lb = bnd;
                
            elseif var_lhs && isnum_rhs && isa_le
                % v <= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.ub = bnd;
                
            elseif var_lhs && isnum_rhs && isa_ge
                % v >= num
                var = obj.find_variable(id_lhs);
                bnd = num_rhs;
                var.lb = bnd;
                
            elseif isnum_lhs && var_rhs && isa_eq
                % num == v
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ub = bnd;
                var.lb = bnd;
                
            elseif isnum_lhs && var_rhs && isa_le
                % num <= v
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.lb = bnd;
                
            elseif isnum_lhs && var_rhs && isa_ge
                % num >= v
                var = obj.find_variable(id_rhs);
                bnd = num_lhs;
                var.ub = bnd;
                
            elseif var_lhs && fnh_rhs && isa_eq
                % v == @
                var = obj.find_variable(id_lhs);
                var.ub = rhs;
                var.lb = rhs;
                
            elseif var_lhs && fnh_rhs && isa_le
                % v <= @
                var = obj.find_variable(id_lhs);
                var.ub = rhs;
                
            elseif var_lhs && fnh_rhs && isa_ge
                % v >= @
                var = obj.find_variable(id_lhs);
                var.lb = rhs;
                
            elseif fnh_lhs && var_rhs && isa_eq
                % @ == v
                var = obj.find_variable(id_rhs);
                var.ub = lhs;
                var.lb = lhs;
                
            elseif fnh_lhs && var_rhs && isa_le
                % @ <= v
                var = obj.find_variable(id_rhs);
                var.lb = lhs;
                
            elseif fnh_lhs && var_rhs && isa_ge
                % @ >= v
                var = obj.find_variable(id_rhs);
                var.ub = lhs;
                
            elseif isa_eq && invariant
                % Transcription invariant equality constraint
                obj.ec_point_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif isa_eq && hard
                % Transcription invariant equality constraint
                obj.ec_hard_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif isa_eq
                % Equality constraint
                obj.ec_eqs{end+1} = canonicalize(ssr).m_lhs;
            
            elseif (isa_le || isa_ge) && invariant
                % Transcription invariant inequality constraint
                obj.iec_point_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif (isa_le || isa_ge) && hard
                % Hard inequality constraint
                obj.iec_hard_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            elseif isa_le || isa_ge
                % Inequality constraint
                obj.iec_eqs{end+1} = canonicalize(ssr).m_lhs;
                
            else
                % error
                error(yop.error.unknown_constraint());
                
            end
        end
        
        function ocp_var = find_variable(obj, id)
            for ocp_var = obj.variables()
                if ocp_var.ast.m_id == id
                    return;
                end
            end
            error(yop.error.failed_to_find_variable(id));
        end
        
        function tp_int = find_special_nodes(obj, expression)
            % Find and add all special nodes of the expression
            
            % Used to find errors in the objective function
            tp_int = {};
            
            % Sort all nodes
            [tsort, N, obj.visited] = ...
                topological_sort(expression, obj.visited);
            
            % Find special nodes, maintain topological order
            for n=1:N
                tn = tsort{n};
                if isa(tn, 'yop.ast_variable')
                    obj.add_variable(tn);
                    
                elseif isa(tn, 'yop.ast_timepoint')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.tp);
                    obj.snodes(end+1) = sn;
                    obj.tps(end+1) = sn;
                    tp_int{end+1} = tn;
                    
                elseif isa(tn, 'yop.ast_int')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.int);
                    obj.snodes(end+1) = sn;
                    obj.ints(end+1) = sn;
                    tp_int{end+1} = tn;
                    
                elseif isa(tn, 'yop.ast_der')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.der);
                    obj.snodes(end+1) = sn;
                    obj.ders(end+1) = sn;
                end
            end
        end
        
    end
    
    %% Canonical form
    methods
        
        function obj = to_canonical(obj)
            yop.progress.ocp_parsing();
            % The system is augmented before box bounds are set in order to
            % use the correct default value for the augemented variables.
            obj.augment_system();
            obj.sort_states();
            obj.set_box_bounds();
            obj.set_objective_fn();
            obj.vectorize_dynamics();
            obj.set_point_con();
            obj.set_path_con();
            obj.set_hard_path_con();
            obj.set_ival_path_con();
            obj.set_special_functions();
            yop.progress.ocp_parsed();
        end
        
        function obj = augment_system(obj)
            
            % Augment system based on control parametrization
            % Step 1: Account for all control inputs
            for uk=obj.controls
                du = uk.ast.m_du;
                while ~isempty(du)
                    obj.add_unique_control(du);
                    du = du.m_du;
                end
            end
            
            % Step 2: Promote integrated controls to states and add
            %         augmenting equations
            keep = [];
            for k=1:length(obj.controls)
                uk = obj.controls(k);
                if ~isempty(uk.ast.m_du)
                    obj.states(end+1) = uk;
                    obj.ode_eqs{end+1} = ode(yop.ast_der(uk.ast)==uk.ast.m_du);
                else
                    keep(end+1) = k;
                end
            end
            obj.controls = obj.controls(keep);

            % Fill in blanks
            if isempty(obj.independent)
                obj.add_independent(yop.independent());
            end
            
            if isempty(obj.independent0)
                obj.add_independent0(yop.independent0());
            end
            
            if isempty(obj.independentf)
                obj.add_independentf(yop.independentf());
            end
            
        end
        
        function obj = set_box_bounds(obj)
            
            % Time0
            if isempty(obj.independent0.ub) && isempty(obj.independent0.lb)
                obj.independent0.ub = yop.defaults.independent0_ub;
                obj.independent0.lb = yop.defaults.independent0_lb;
            end
                
            if isempty(obj.independent0.ub)
                obj.independent0.ub = inf;
            end
                
            if isempty(obj.independent0.lb)
                obj.independent0.lb = -inf;
            end
            
            % Timef
            if isempty(obj.independentf.ub) && isempty(obj.independentf.lb)
                obj.independentf.ub = yop.defaults.independentf_ub;
                obj.independentf.lb = yop.defaults.independentf_lb; 
            end
            
            if isempty(obj.independentf.ub)
                obj.independentf.ub = inf;
            end
            
            if isempty(obj.independentf.lb)
                % User should introduce t0 <= tf
                warning(['[Yop] Final time is unbounded from below. ', ...
                    'If this is intentional, consider introducing ''t0 <= tf''. ', ...
                    'If this was unintentional, set a lower bound for tf ' ...
                    '(''tf >= value'') and consider introduction the above mentioned ' ...
                    'constraint if t0 and tf can overlap.']);
                obj.independentf.lb = -inf; 
            end
            
            % Time
            % Enables constraints such as 't > 0, t < 10'
            if ~isempty(obj.independent.lb)
                t_min = obj.independent.lb;
                obj.independent0.lb = max(obj.independent0.lb, t_min);
                obj.independent0.ub = max(obj.independent0.ub, t_min);
                obj.independentf.lb = max(obj.independentf.lb, t_min);
                obj.independentf.ub = max(obj.independentf.ub, t_min);
            end
            
            if ~isempty(obj.independent.lb)
                t_max = obj.independent.ub;
                obj.independent0.lb = min(obj.independent0.lb, t_max);
                obj.independent0.ub = min(obj.independent0.ub, t_max);
                obj.independentf.lb = min(obj.independentf.lb, t_max);
                obj.independentf.ub = min(obj.independentf.ub, t_max);
            end
            
            % State
            for x=obj.states
                if isempty(x.ub)
                    x.ub = yop.defaults.state_ub;
                end
                if isempty(x.lb)
                    x.lb = yop.defaults.state_lb;
                end
                if isempty(x.ub0)
                    x.ub0 = x.ub;
                end
                if isempty(x.lb0)
                    x.lb0 = x.lb; 
                end
                if isempty(x.ubf)
                    x.ubf = x.ub;
                end
                if isempty(x.lbf)
                    x.lbf = x.lb;
                end
            end
            
            % Algebraics 
            for z=obj.algebraics
                if isempty(z.ub)
                    z.ub = yop.defaults.algebraic_ub;
                end
                if isempty(z.lb)
                    z.lb = yop.defaults.algebraic_lb;
                end
            end
            
            % Controls
            for u=obj.controls
                if isempty(u.ub)
                    u.ub = yop.defaults.control_ub;
                end
                if isempty(u.lb)
                    u.lb = yop.defaults.control_lb;
                end
                if isempty(u.ub0)
                    u.ub0 = u.ub;
                end
                if isempty(u.lb0)
                    u.lb0 = u.lb; 
                end
                if isempty(u.ubf)
                    u.ubf = u.ub;
                end
                if isempty(u.lbf)
                    u.lbf = u.lb;
                end
            end
            
            % Parameters
            for p=obj.parameters
                if isempty(p.ub)
                    p.ub = yop.defaults.parameter_ub;
                end
                if isempty(p.lb)
                    p.lb = yop.defaults.parameter_lb;
                end
            end
        end
        
        function obj = set_objective_fn(obj)
            t0 = mx_vec(obj.independent0);
            tf = mx_vec(obj.independentf);
            pp = mx_vec(obj.parameters);
            tp = mx_vec(obj.tps);
            ii = mx_vec(obj.ints);
            args = {t0,tf,pp,tp,ii};
            J = casadi.Function('J', args, {value(obj.objective.ast)});
            
            % The nlp variables are scaled, so they must be descaled before
            % the expression can be evaluate.
            t0s = t0 .* obj.W_t0 - obj.OS_t0;
            tfs = tf .* obj.W_tf - obj.OS_tf;
            pps = pp .* obj.W_p  - obj.OS_p;
            Js  = J(t0s, tfs, pps, tp, ii);
            obj.objective.fn = casadi.Function('Js', args, {Js});
        end 
        
        function obj = set_special_functions(obj)
            for sn = obj.snodes
                sn.fn = obj.dsfn(value(sn.ast.m_expr));
            end
        end
        
        function obj = vectorize_dynamics(obj)
            % Vectorize the equation
            n_ode = length(obj.ode_eqs);
            tmp_lhs = cell(n_ode, 1);
            tmp_rhs = cell(n_ode, 1);
            for k=1:length(obj.ode_eqs)
                tmp_lhs{k} = obj.ode_eqs{k}.m_lhs;
                tmp_rhs{k} = obj.ode_eqs{k}.m_rhs;
            end
            ode_lhs = vertcat(tmp_lhs{:});
            
            % Test if all states are bound to an ode
            [~, ode_ids] = Type(ode_lhs);
            [ode_ids, idx] = sort(ode_ids);
            x_ids = obj.get_state_ids();
            if ~isequal(x_ids, ode_ids)
                state_ast = {};
                for id = setdiff(ode_ids, x_ids)
                    state_ast{end+1} = obj.find_variable(id);
                end
                error(yop.error.missing_state_derivative(state_ast));
            end
            
            % Change order of equations so that state vector and ode
            % equation order match
            obj.ode.m_lhs = ode_lhs(idx);
            obj.ode.m_rhs = vertcat(tmp_rhs{idx});
            
            % Algebraic equation
            nz = length(obj.alg_eqs);
            alg_rhs = cell(nz,1);
            for k=1:nz
                z_k = obj.alg_eqs{k};
                alg_rhs{k} = z_k.m_rhs - z_k.m_lhs;
            end
            obj.alg.m_lhs = zeros(nz, 1);
            obj.alg.m_rhs = vertcat(alg_rhs{:});
            
            % Compute symbolic functions of the dynamics
            obj.set_dynamics_fn();
        end
        
        function obj = set_dynamics_fn(obj)
            ode_expr = value(obj.ode.m_rhs);
            alg_expr = value(obj.alg.m_rhs);
            f = casadi.Function('f', obj.mx_args(), {ode_expr});
            
            % Descale input variables, evaluate ode rhs, scale derivative
            dargs = obj.mx_dargs();
            fs = f(dargs{:}).*(1./obj.W_x);
            obj.ode.fn = casadi.Function('ode', obj.mx_args(), {fs});
            obj.alg.fn = obj.dsfn(alg_expr, 'alg');
        end
        
        function obj = set_path_con(obj)                       
            expr = value(vertcat(obj.ec_eqs{:}, obj.iec_eqs{:}));
            obj.path.fn = obj.dsfn(expr(:), 'eq');
            n_eq  = length(obj.ec_eqs);
            n_ieq = length(obj.iec_eqs);
            obj.path.ub = zeros(n_eq+n_ieq, 1);
            obj.path.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function obj = set_hard_path_con(obj)
            expr = value(vertcat(obj.ec_hard_eqs{:}, obj.iec_hard_eqs{:}));
            obj.path_hard.fn = obj.dsfn(expr(:), 'heq');
            n_eq  = length(obj.ec_hard_eqs);
            n_ieq = length(obj.iec_hard_eqs);
            obj.path_hard.ub = zeros(n_eq+n_ieq, 1);
            obj.path_hard.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function obj = set_point_con(obj)
            expr = value(vertcat(obj.ec_point_eqs{:}, obj.iec_point_eqs{:}));
            obj.point.fn = obj.dsfn(expr(:), 'peq');
            n_eq = length(obj.ec_point_eqs);
            n_ieq = length(obj.iec_point_eqs);
            obj.point.ub = zeros(n_eq+n_ieq, 1);
            obj.point.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function obj = set_ival_path_con(obj)
            data = yop.ival.empty(1,0);
            for k=1:length(obj.ec_ival_eqs)
                ik = obj.ec_ival_eqs{k};
                [t0, tf] = get_ival(ik);
                data(k).t0 = t0;
                data(k).tf = tf;
                data(k).ast = ik;
                data(k).ub = 0;
                data(k).lb = 0;
                expr = value(ik);
                data(k).fn = obj.dsfn(expr(:), 'iveq');
            end
            
            for k=1:length(obj.iec_ival_eqs)
                ik = obj.iec_ival_eqs{k};
                [t0, tf] = get_ival(ik);
                data(k).t0 = t0;
                data(k).tf = tf;
                data(k).ast = ik;
                data(k).ub = 0;
                data(k).lb = -inf;
                expr = value(ik);
                data(k).fn = obj.dsfn(expr(:), 'ivieq');
            end
            
            obj.path_ival = data;
        end
        
        function args = mx_vars(obj)
            args = { ...
                mx_vec(obj.independent0), ...
                mx_vec(obj.independentf), ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters) ...
                };
        end
        
        function args = mx_args(obj)
            args = { ...
                mx_vec(obj.independent0), ...
                mx_vec(obj.independentf), ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters), ...
                mx_vec(obj.tps), ...
                mx_vec(obj.ints), ...
                mx_vec(obj.ders) ...
                };
        end
        
        function fn = dsfn(obj, expr, fname)
            % descaled function
            if nargin == 2
                fname = 'f';
            end
            args  = obj.mx_args();
            dargs = obj.mx_dargs();
            tmp = casadi.Function(fname, args, {expr});
            fn  = casadi.Function(fname, args, {tmp(dargs{:})});
        end
        
        function args = mx_dargs(obj)
            % Descale
            t0 = obj.descale_t0(mx_vec(obj.independent0));
            tf = obj.descale_tf(mx_vec(obj.independentf));
            tt = obj.descale_t (mx_vec(obj.independent));
            xx = obj.descale_x (mx_vec(obj.states));
            zz = obj.descale_z (mx_vec(obj.algebraics));
            uu = obj.descale_u (mx_vec(obj.controls));
            pp = obj.descale_p (mx_vec(obj.parameters));
            tp = mx_vec(obj.tps);
            ii = mx_vec(obj.ints);
            dd = mx_vec(obj.ders);
            args = {t0, tf, tt, xx, zz, uu, pp, tp, ii, dd};
        end
        
        function vs = scale_t0(obj, v)
            % scale t0
            vs = (v + obj.OS_t0).*(1./obj.W_t0);
        end
        
        function vs = scale_tf(obj, v)
            % scale tf
            vs = (v + obj.OS_tf).*(1./obj.W_tf);
        end
        
        function vs = scale_t(obj, v)
            % scale t
            vs = (v + obj.OS_t).*(1./obj.W_t);
        end
        
        function vs = scale_x(obj, v)
            % scale x
            vs = (v + obj.OS_x).*(1./obj.W_x);
        end
        
        function vs = scale_z(obj, v)
            % scale z
            vs = (v + obj.OS_z).*(1./obj.W_z);
        end
        
        function vs = scale_u(obj, v)
            % scale u
            vs = (v + obj.OS_u).*(1./obj.W_u);
        end
        
        function vs = scale_p(obj, v)
            % scale p
            vs = (v + obj.OS_p).*(1./obj.W_p);
        end
        
        function v = descale_t0(obj, vs)
            % descale t0
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_t0 - obj.OS_t0;
            end
        end
        
        function v = descale_tf(obj, vs)
            % descale tf
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_tf - obj.OS_tf;
            end
        end
        
        function v = descale_t(obj, vs)
            % descale t
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_t - obj.OS_t;
            end
        end
        
        function v = descale_x(obj, vs)
            % descale x
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_x - obj.OS_x;
            end
        end
        
        function v = descale_z(obj, vs)
            % descale z
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_z - obj.OS_z;
            end
        end
        
        function v = descale_u(obj, vs)
            % descale u
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_u - obj.OS_u;
            end
        end
        
        function v = descale_p(obj, vs)
            % descale p
            if isempty(vs)
                v = vs;
            else
                v = vs.*obj.W_p - obj.OS_p;
            end
        end
        
        function W = W_t0(obj)
            W = obj.independent0.weight;
        end
        
        function W = W_tf(obj)
            W = obj.independentf.weight;
        end
        
        function W = W_t(obj)
            W = obj.independent.weight;
        end
        
        function W = W_x(obj)
            W = obj.states.weight;
        end
        
        function W = W_z(obj)
            W = obj.algebraics.weight;
        end
        
        function W = W_u(obj)
            W = obj.controls.weight;
        end
        
        function W = W_p(obj)
            W = obj.parameters.weight;
        end
        
        function OS = OS_t0(obj)
            OS = obj.independent0.offset;
        end
        
        function OS = OS_tf(obj)
            OS = obj.independentf.offset;
        end
        
        function OS = OS_t(obj)
            OS = obj.independent.offset;
        end
        
        function OS = OS_x(obj)
            OS = obj.states.offset;
        end
        
        function OS = OS_z(obj)
            OS = obj.algebraics.offset;
        end
        
        function OS = OS_u(obj)
            OS = obj.controls.offset;
        end
        
        function OS = OS_p(obj)
            OS = obj.parameters.offset;
        end
        
        function IDs = ids(obj)
            IDs = [obj.independent0.ids, obj.independentf.ids, ...
                obj.independent.ids, obj.states.ids, obj.algebraics.ids, ...
                obj.controls.ids, obj.parameters.ids];
        end
        
        function ids = get_state_ids(obj)
            nx = length(obj.states);
            ids = zeros(nx,1);
            for k=1:nx
                ids(k) = obj.states(k).ast.m_id;
            end
        end
        
        function ids = sort_states(obj)
            % Change state order so that they come in id order
            [ids, idx] = sort(obj.get_state_ids());
            obj.states = obj.states(idx);
        end
        
        
        
        function v = variables(obj)
            v = [obj.independent, obj.independent0, obj.independentf, ...
                obj.states, obj.algebraics, obj.controls, obj.parameters];
        end
        
        function remove_state_der(obj, id)
            % Remove the state derivative node if all elements covered by
            % the derivative are states.
            for k=1:length(obj.ders)
                if all(isa_der(obj.ders(k).ast)) && all(obj.ders(k).ast.m_der == id)
                    if all(Type(obj.ders(k).ast) == yop.var_type.state)
                        to_remove = obj.ders(k);
                        obj.ders = [obj.ders(1:k-1), obj.ders(k+1:end)]; 
                        % Also need to remove it from special nodes vector
                        for n=1:length(obj.snodes)
                            if obj.snodes(n) == to_remove
                                obj.snodes = [obj.snodes(1:n-1), ...
                                    obj.snodes(n+1:end)];
                                return;
                            end
                        end
                    end
                end
            end
        end
        
    end
       
    %% Helper functions / Misc
    methods
        function add_variable(obj, v)
            switch class(v)
                case 'yop.ast_independent'
                    obj.add_independent(v);
                case 'yop.ast_independent_initial'
                    obj.add_independent0(v);
                case 'yop.ast_independent_final'
                    obj.add_independentf(v);
                case 'yop.ast_state'
                    obj.add_state(v);
                case 'yop.ast_algebraic'
                    obj.add_algebraic(v);
                case 'yop.ast_control'
                    obj.add_control(v);
                case 'yop.ast_parameter'
                    obj.add_parameter(v);
            end
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ocp_var(t);
            else
                yop.error.multiple_independent_variables();
            end
        end
        
        function obj = add_independent0(obj, t)
            if isempty(obj.independent0)
                obj.independent0 = yop.ocp_var(t);
            else
                error(yop.error.multiple_independent_initial());
            end
        end
        
        function obj = add_independentf(obj, t)
            if isempty(obj.independentf)
                obj.independentf = yop.ocp_var(t);
            else
                error(yop.error.multiple_independent_final());
            end
        end
        
        function obj = add_state(obj, x)
            obj.states(end+1) = yop.ocp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            obj.algebraics(end+1) = yop.ocp_var(z);
        end
        
        function obj = add_control(obj, u)
            obj.controls(end+1) = yop.ocp_var(u);
        end
        
        function obj = add_unique_control(obj, u)
            for uk = obj.controls
                if uk.ast.m_id == u.m_id
                    return;
                end
            end
            obj.add_control(u);
        end
        
        function obj = add_parameter(obj, p)
            obj.parameters(end+1) = yop.ocp_var(p);
        end
    end
    
    %% NLP Interface
    methods
        
        function n = n_x(obj)
            n = length(obj.states);
        end
        
        function n = n_z(obj)
            n = length(obj.algebraics);
        end
        
        function n = n_u(obj)
            n = length(obj.controls);
        end
        
        function n = n_p(obj)
            n = length(obj.parameters);
        end
        
        function n = n_tp(obj)
            n = n_elem(obj.tps);
        end
        
        function n = n_int(obj)
            n = n_elem(obj.ints);
        end
        
        function n = n_der(obj)
            n = n_elem(obj.ders);
        end
        
        function bool = has_initial_guess(obj)
            bool = ~isempty(obj.guess);
        end
        
        function bool = has_path(obj)
            bool = ~isempty(obj.path.ub);
        end
        
        function bool = has_hard_path(obj)
            bool = ~isempty(obj.path_hard.ub);
        end
        
        function bd = t0_ub(obj, t)
            if isa(obj.independent0.ub, 'function_handle')
                bd = obj.independent0.ub(t);
            else
                bd = obj.independent0.ub;
            end
            bd = obj.scale_t0(bd);
        end
        
        function bd = t0_lb(obj, t)
            if isa(obj.independent0.lb, 'function_handle')
                bd = obj.independent0.lb(t);
            else
                bd = obj.independent0.lb;
            end
            bd = obj.scale_t0(bd);
        end
        
        function bd = tf_ub(obj, t)
            if isa(obj.independentf.ub, 'function_handle')
                bd = obj.independentf.ub(t);
            else
                bd = obj.independentf.ub;
            end
            bd = obj.scale_tf(bd);
        end
        
        function bd = tf_lb(obj, t)
            if isa(obj.independentf.lb, 'function_handle')
                bd = obj.independentf.lb(t);
            else
                bd = obj.independentf.lb;
            end
            bd = obj.scale_tf(bd);
        end
        
        function bd = x0_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ub0, 'function_handle')
                    bd(end+1) = v.ub0(t);
                else
                    bd(end+1) = v.ub0;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = x0_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lb0, 'function_handle')
                    bd(end+1) = v.lb0(t);
                else
                    bd(end+1) = v.lb0;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = x_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ub, 'function_handle')
                    bd(end+1) = v.ub(t);
                else
                    bd(end+1) = v.ub;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = x_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lb, 'function_handle')
                    bd(end+1) = v.lb(t);
                else
                    bd(end+1) = v.lb;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = xf_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ubf, 'function_handle')
                    bd(end+1) = v.ubf(t);
                else
                    bd(end+1) = v.ubf;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = xf_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lbf, 'function_handle')
                    bd(end+1) = v.lbf(t);
                else
                    bd(end+1) = v.lbf;
                end
            end
            bd = obj.scale_x(bd(:));
        end
        
        function bd = z_ub(obj, t)
            bd = [];
            if ~isempty(obj.algebraics)
                if isa(obj.algebraics.ub, 'function_handle')
                    bd = obj.algebraics.ub(t);
                else
                    bd = obj.algebraics.ub;
                end
            end
            bd = obj.scale_z(bd(:));
        end
        
        function bd = z_lb(obj, t)
            bd = [];
            if ~isempty(obj.algebraics)
                if isa(obj.algebraics.lb, 'function_handle')
                    bd = obj.algebraics.lb(t);
                else
                    bd = obj.algebraics.lb;
                end
            end
            bd = obj.scale_z(bd(:));
        end
        
        function bd = u0_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ub0, 'function_handle')
                    bd(end+1) = v.ub0(t);
                else
                    bd(end+1) = v.ub0;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = u0_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lb0, 'function_handle')
                    bd(end+1) = v.lb0(t);
                else
                    bd(end+1) = v.lb0;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = u_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ub, 'function_handle')
                    bd(end+1) = v.ub(t);
                else
                    bd(end+1) = v.ub;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = u_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lb, 'function_handle')
                    bd(end+1) = v.lb(t);
                else
                    bd(end+1) = v.lb;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = uf_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ubf, 'function_handle')
                    bd(end+1) = v.ubf(t);
                else
                    bd(end+1) = v.ubf;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = uf_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lbf, 'function_handle')
                    bd(end+1) = v.lbf(t);
                else
                    bd(end+1) = v.lbf;
                end
            end
            bd = obj.scale_u(bd(:));
        end
        
        function bd = p_ub(obj, t)
            bd = [];
            if ~isempty(obj.parameters)
                if isa(obj.parameters.ub, 'function_handle')
                    bd = obj.parameters.ub(t);
                else
                    bd = obj.parameters.ub;
                end
            end
            bd = obj.scale_p(bd(:));
        end
        
        function bd = p_lb(obj, t)
            bd = [];
            if ~isempty(obj.parameters)
                if isa(obj.parameters.lb, 'function_handle')
                    bd = obj.parameters.lb(t);
                else
                    bd = obj.parameters.lb;
                end
            end
            bd = obj.scale_p(bd(:));
        end
        
        function [t0, tf] = get_horizon(obj)
            t0_ub = obj.independent0.ub;
            t0_lb = obj.independent0.lb;
            tf_ub = obj.independentf.ub;
            tf_lb = obj.independentf.lb;
            
            fixed = t0_lb == t0_ub && tf_lb == tf_ub && ...
                all(~isinf([t0_lb, t0_ub, tf_lb, tf_ub]));
            
            if fixed
                t0 = t0_lb;
                tf = tf_lb;
            else
                % Normalized value
                t0 = 0;
                tf = 1;
            end
        end
        
        function [t00, tf0, t0, x0, z0, u0, p0] = initial_guess(obj)
            t0  = obj.guess.value(obj.independent.ast)'; % Do not scale, used for interpolation only
            t00 = obj.scale_t0(t0(1));
            tf0 = obj.scale_tf(t0(end));
            x0  = obj.state_guess(t0);
            z0  = obj.algebraic_guess(t0);
            u0  = obj.control_guess(t0);
            p0  = obj.scale_p(obj.guess.p(:));
        end
        
        function x0 = state_guess(obj, t0)
            x0 = [];
            W = obj.W_x;
            OS = obj.OS_x;
            for x = obj.states
                xk = obj.guess.value(x.ast);
                if isempty(xk)
                    xk = ones(length(t0), 1);
                end
                x0 = [x0, xk(:)];
            end
            x0 = (x0 + obj.OS_x').*(1./obj.W_x');
        end
           
        function z0 = algebraic_guess(obj, t0)
            z0 = [];
            for z = obj.algebraics
                zk = obj.guess.value(z.ast);
                if isempty(zk)
                    zk = ones(length(t0), 1);
                end
                z0 = [z0, zk(:)];
            end
            z0 = (z0 + obj.OS_z').*(1./obj.W_z');
        end
        
        function u0 = control_guess(obj, t0)
            u0 = [];
            for u = obj.controls
                uk = obj.guess.value(u.ast);
                if isempty(uk)
                    uk = ones(length(t0), 1);
                end
                u0 = [u0, uk(:)];
            end
            u0 = (u0 + obj.OS_u').*(1./obj.W_u');
        end
        
    end
    
    
    %% Multi-phase - Overloading for multiple phase problems
    methods
        function sum = plus(lhs, rhs)
            sum = yop.ocp('');
            sum.m_value = lhs.m_value + rhs.m_value;
            sum.m_phases = [lhs.m_phases, rhs];
            % parent is only valid for sequential problems
            sum.m_parent = [lhs.m_parent, lhs.m_phases(end)];
        end
    end
    
    %% Static
    methods (Static)
        function ID = get_id()
            persistent cur
            if isempty(cur)
                cur = 0;
            end
            cur = cur + 1;
            ID = cur;
        end
    end
end