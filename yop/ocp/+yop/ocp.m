classdef ocp < handle
    properties
        % Misc
        name
        
        % objective
        objective
        
        % Variables
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        differetial_eqs
        constraints
        
        % Vectorized form
        Px % State perturbation function
        t0
        tf
        t
        x
        z
        u
        p
        ode
        
        % Vectorized bounds
        t0_ub
        t0_lb
        tf_ub
        tf_lb
                
        x0_ub
        x0_lb
        x_ub
        x_lb
        xf_ub
        xf_lb
        
        z_ub
        z_lb
        
        u0_ub
        u0_lb
        u_ub
        u_lb
        uf_ub
        uf_lb
        
        p_ub
        p_lb
        
        % Parsed constraints
        eq = {}; % Temporary cell
        ieq = {}; % Temporary cell
    end
    
    %% Formulate ocp
    methods 
        
        function obj = ocp(name)
            if nargin == 1
                obj.name = name;
            else
                obj.name = '';
            end
        end
        
        function obj = min(obj, objective)
            obj.objective = objective;
        end
        
        function obj = max(obj, objective)
            % Maybe its better to note that it is a max problem, so that
            % the present function can present it as a maximization problem
            obj.objective = -objective;
        end
        
        function obj = st(obj, varargin)
            obj.constraints = varargin;
        end
        
    end
    
    %% Analyze user input, and build canonical form
    methods 

        function obj = build(obj)
            % 1) Find all variables from the relations and expressions that
            %    make up the problem
            obj.parse_variables();
            
            % 2) Convert to single relation form
            srf = yop.to_srf(obj.constraints);
            
            % 3) Convert to homogeneous srf
            hsrf = yop.to_hsrf(srf.get_relations());
            
            % 4) Convert to value-numeric form
            vnf = yop.to_vnf(hsrf);
            
            % 5) Convert to distinct timepoint form (final form for boxcon)
            dtp = yop.to_dtp(vnf);
            obj.parse_box_constraints(dtp);
            obj.set_box_bounds();
            
            % 6) Sort the nonbox constraints: ode, inequality, equality
            nbc = yop.sort_nonbox(dtp.get_pathcon());
            obj.differetial_eqs = nbc.odes;
            
            % 7) Convert to transcription invariant/variant form
            tri = yop.tr_invar(nbc); 
            obj.eq = {tri.eq_inv{:}, tri.eq_var{:}};
            obj.ieq = {tri.ieq_inv{:}, tri.ieq_var{:}};
            
            % 8) vectorize ocp
            obj.vectorize();
            
            % 9) Error checking
            
            % 10) OCP built
                        
        end
        
        function obj = parse_variables(obj)
            vars = yop.get_vars({obj.objective, obj.constraints{:}});
            for k=1:length(vars)
                vk = vars{k};
                switch class(vk)
                    case 'yop.ast_independent'
                        obj.add_independent(vk);
                    case 'yop.ast_independent_initial'
                        obj.add_independent_initial(vk);
                    case 'yop.ast_independent_final'
                        obj.add_independent_final(vk);
                    case 'yop.ast_state'
                        obj.add_state(vk);
                    case 'yop.ast_algebraic'
                        obj.add_algebraic(vk);
                    case 'yop.ast_control'
                        obj.add_control(vk);
                    case 'yop.ast_parameter'
                        obj.add_parameter(vk);
                    otherwise
                        error('[Yop] Error: Unknown variable type');
                end
            end
            
            % If we do not have some of the variables, they are set to
            % something in order for the transcription methods to query
            % them, and also to be able to set bound on independent
            % variables.
            if isempty(obj.independent)
                obj.add_independent(yop.ast_independent('t'));
            end
            
            if isempty(obj.independent_initial)
                obj.add_independent_initial( ...
                    yop.ast_independent_initial('t0'));
            end
            
            if isempty(obj.independent_final)
                obj.add_independent_final( ...
                    yop.ast_independent_final('tf'));
            end
            
            if isempty(obj.states)
                obj.add_state(yop.ast_state('x', 0, 0));
            end
            
            if isempty(obj.algebraics)
                obj.add_algebraic(yop.ast_algebraic('z', 0, 0));
            end
            
            if isempty(obj.controls)
                obj.add_control(yop.ast_control('u', 0, 0));
            end
            
            if isempty(obj.parameters)
                obj.add_parameter(yop.ast_parameter('p', 0, 0));
            end
            
        end
        
        function obj = parse_box_constraints(obj, dtp)
            % Box constraints: var-num
            for k=1:length(dtp.vn_t)
                bc = dtp.vn_t{k};
                var_expr = bc.lhs;
                bnd = yop.prop_num(bc.rhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ub(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lb(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: var-num, t==t0
            for k=1:length(dtp.vn_t0)
                bc = dtp.vn_t0{k};
                var_expr = bc.lhs;
                bnd = yop.prop_num(bc.rhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ub0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lb0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: var-num, t==tf
            for k=1:length(dtp.vn_tf)
                bc = dtp.vn_tf{k};
                var_expr = bc.lhs;
                bnd = yop.prop_num(bc.rhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ubf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lbf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ubf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lbf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var
            for k=1:length(dtp.nv_t)
                bc = dtp.nv_t{k};
                var_expr = bc.rhs;
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ub(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lb(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var, t==t0
            for k=1:length(dtp.nv_t0)
                bc = dtp.nv_t0{k};
                var_expr = bc.rhs;
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ub0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lb0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub0(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var, t==tf
            for k=1:length(dtp.nv_tf)
                bc = dtp.nv_tf{k};
                var_expr = bc.rhs;
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ubf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                        var.lbf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lbf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ubf(re.idx_var) = yop.get_subexpr(bnd, re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
        end
        
        
        function obj = set_box_bounds(obj)
            obj.set_box_bnd('independent', 'lb', 'ub', ...
                yop.defaults().independent_lb, ...
                yop.defaults().independent_ub);
            
            obj.set_box_bnd('independent_initial', 'lb', 'ub', ...
                yop.defaults().independent_lb0, ...
                yop.defaults().independent_ub0);
            
            obj.set_box_bnd('independent_final', 'lb', 'ub', ...
                yop.defaults().independent_lbf, ...
                yop.defaults().independent_ubf);
            
            obj.set_box_bnd('states', 'lb', 'ub', ...
                yop.defaults().state_lb, yop.defaults().state_ub);
            
            obj.set_box_bnd('states', 'lb0', 'ub0', ...
                yop.defaults().state_lb0, yop.defaults().state_ub0);
            
            obj.set_box_bnd('states', 'lbf', 'ubf', ...
                yop.defaults().state_lbf, yop.defaults().state_ubf);
            
            obj.set_box_bnd('algebraics', 'lb', 'ub', ...
                yop.defaults().algebraic_lb, yop.defaults().algebraic_ub);
            
            obj.set_box_bnd('controls', 'lb', 'ub', ...
                yop.defaults().control_lb, yop.defaults().control_ub);
            
            obj.set_box_bnd('controls', 'lb0', 'ub0', ...
                yop.defaults().control_lb0, yop.defaults().control_ub0);
            
            obj.set_box_bnd('controls', 'lbf', 'ubf', ...
                yop.defaults().control_lbf, yop.defaults().control_ubf);
            
            obj.set_box_bnd('parameters', 'lb', 'ub', ...
                yop.defaults().parameter_lb, yop.defaults().parameter_ub);
        end
        
        function obj = set_box_bnd(obj, var_field, lb_field, ...
                ub_field, lb_def, ub_def)
            
            for v=obj.(var_field)
                % Timed box constraint: t=t0 or t=tf, follow a hierarchical
                % approach when it comes to setting the. Beginning with the
                % most important the logic is:
                %   1) A value set by the user has highest precedence.
                %   2) If a timed value is not set, but one is set for the
                %      horizon t=(t0, tf) then that value is used.
                %   3) If that too is unset, then the default value is
                %      used.
                
                ub = v.(ub_field); % The full bound vector
                not_set = isnan(ub); % The elements not set
                
                bd = v.ub(not_set); % The ub for t = (t0, tf)
                bd(isnan(bd)) = ub_def; % Default for unsets 
                
                ub(not_set) = bd; % Timed bound
                v.(ub_field) = ub; % Set variable timed bound
                
                % Same procedure for lb
                lb = v.(lb_field);
                not_set = isnan(lb);
                bd = v.lb(not_set);
                bd(isnan(bd)) = lb_def;
                lb(not_set) = bd;
                v.(lb_field) = lb;
            end
        end
        
        function var = find_variable(obj, id)
            for v = obj.variables
                if v.var.id == id
                    var = v;
                    return;
                end
            end
            error('[Yop] Error: ID not found');
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'independent variable']);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'initial bound for the independent variable']);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'terminal bound for the independent variable']);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ocp_var(x);
            else
                obj.states(end+1) = yop.ocp_var(x);
            end
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_var(z);
            else
                obj.algebraics(end+1) = yop.ocp_var(z);
            end
        end
        
        function obj = add_algebraic_eq(obj, eq)
            if isempty(obj.alg_eqs)
                obj.alg_eqs = yop.ocp_var(eq);
            else
                obj.alg_eqs(end+1) = yop.ocp_var(eq);
            end
        end
        
        function obj = add_control(obj, u)
            if isempty(obj.controls)
                obj.controls = yop.ocp_var(u);
            else
                obj.controls(end+1) = yop.ocp_var(u);
            end
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_var(p);
            else
                obj.parameters(end+1) = yop.ocp_var(p);
            end
        end
        
        function vars = variables(obj)
            vars = [obj.independent(:).', ...
                obj.independent_initial(:).', ...
                obj.independent_final(:).', ...
                obj.states(:).', ...
                obj.algebraics(:).', ...
                obj.controls(:).', ...
                obj.parameters(:).'];
        end
        
        function obj = vectorize(obj)
            obj.vectorize_independent();
            obj.vectorize_state();
            obj.vectorize_algebraic();
            obj.vectorize_control();
            obj.vectorize_parameters();
        end
        
        function obj = vectorize_independent(obj)
            obj.t = obj.independent_initial.mx;
            obj.t0 = obj.independent_initial.mx;
            obj.tf = obj.independent_final.mx;
            
            obj.t0_lb = obj.independent_initial.lb;
            obj.t0_ub = obj.independent_initial.ub;
            obj.tf_lb = obj.independent_final.lb;
            obj.tf_ub = obj.independent_final.ub;
        end
        
        function obj = vectorize_state(obj)
            %  VECTORIZE_STATE - Vectorize state and ode
            %  This function is the reason why the ocp needs to be
            %  vectorized. The odes governing the states may have have the
            %  variable vector (lhs of 'der(x) == f(t,x,z,u,p)') in any 
            %  order (x can be permuted or lack some of the sates on the 
            %  lhs). This becomes a problem in integration, where the
            %  righthandside is evaluated and added to the current step in
            %  order to compute the next (x_next = x_cur + h*f(.)). Unless
            %  the order match, the integration will not be correct. For
            %  that reason it is necessary to permute the state vector, x,
            %  so that the elements come in the same order as their
            %  derivatives.
            
            %  1) Vectorize ode expression
            ode_var = obj.differetial_eqs(1).var;
            ode_expr = obj.differetial_eqs(1).expr;
            for k=2:length(obj.differetial_eqs)
                ode_var = [ode_var(:); obj.differetial_eqs(k).var(:)];
                ode_expr = [ode_expr(:); obj.differetial_eqs(k).expr(:)];
            end
            
            assert(all(isa_variable(ode_var)), '[Yop] Unexpected error.');
            
            %  2) Perform a reaching elements analysis in order to
            %  determine what variables, and what exact elements reaches
            %  the der(...) expression. 
            x_cell = arrayfun(@(e) e.var, obj.states, ...
                'UniformOutput', false);
            [re, nr] = yop.reaching_elems(ode_var, x_cell);
            
            %  3) Three cases: 
            %        i) All elements of the variable reach (nothing to do).
            %       ii) Some of the elements reach (add not reaching).
            %      iii) None of the elements reach (add all).
            % case i && ii
            for k=1:length(re)
                if numel(re(k).enum) ~= numel(re(k).reaching)
                    % Number of elements and the number reaching does not
                    % match. Pick out the ones not reaching and add them to
                    % the variable vector and set their derivative to zero.
                    [~, idx] = setdiff(re(k).enum, re(k).reaching);
                    vk = re(k).var(idx);
                    ode_var = [ode_var(:); vk(:)];
                    ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                end
            end
            %
            % case iii) do not reach at all, so they are added in full.
            for k=1:length(nr)
                vk = nr(k).var;
                ode_var = [ode_var(:); vk(:)];
                ode_expr = [ode_expr(:); zeros(size(vk(:)))];
            end
            
            % 4) Create a perturbation function, that takes the
            % concatenated state variables and perturbes them to mathch the
            % vectorized ode. It is necessary to complicate this process a
            % bit, since the ast that make up the lhs may contain variables
            % that never reach the end, there is a risk that these are free
            % variables on the casadi function, so to avoid this, the
            % perturbation function is created in two steps. First,
            % deriving an expression for the lhs, given all arguments as
            % input. And in the second step use the first function to
            % derive a function that only takes the concatenated states as
            % input.
            obj.set_mx(); % use of casadi for convenience, other ways poss.
            [~, ~, t, x, u, p] = obj.mx_vecs();
            x_fn = casadi.Function('x', {t,x,u,p}, {evaluate(ode_var)});
            P_expr = x_fn(t, x, u, p); % perturbation expression
            obj.Px = casadi.Function('P', {x}, {P_expr}); % Perturbation fn
            obj.x = obj.Px(x);
            
            % 5) Set ode_rhs
            ode_fn = casadi.Function('x', {t,obj.x,u,p}, {evaluate(ode_expr)});
            obj.ode = ode_fn;
            
            % 6) Set state bounds
            x0_lb=[]; x0_ub=[]; x_lb=[]; x_ub=[]; xf_lb=[]; xf_ub=[];
            for k=1:length(obj.states)
                x_lb = [x_lb(:); obj.states(k).lb(:)];
                x_ub = [x_ub(:); obj.states(k).ub(:)];
                x0_lb = [x0_lb(:); obj.states(k).lb0(:)];
                x0_ub = [x0_ub(:); obj.states(k).ub0(:)];
                xf_lb = [xf_lb(:); obj.states(k).lbf(:)];
                xf_ub = [xf_ub(:); obj.states(k).ubf(:)];
            end
            obj.x0_lb = full(obj.Px(x0_lb));
            obj.x0_ub = full(obj.Px(x0_ub));
            obj.x_lb  = full(obj.Px(x_lb));
            obj.x_ub  = full(obj.Px(x_ub));
            obj.xf_lb = full(obj.Px(xf_lb));
            obj.xf_ub = full(obj.Px(xf_ub));
        end
        
        function obj = vectorize_algebraic(obj)
            obj.z = obj.algebraics.mx_vec();
            z_lb=[]; z_ub=[];
            for k=1:length(obj.algebraics)
                z_lb = [z_lb(:); obj.algebraics(k).lb(:)];
                z_ub = [z_ub(:); obj.algebraics(k).ub(:)];
            end
            obj.z_lb  = z_lb;
            obj.z_ub  = z_ub;
        end
        
        function obj = vectorize_control(obj)
            obj.u = obj.controls.mx_vec();
            u0_lb=[]; u0_ub=[]; u_lb=[]; u_ub=[]; uf_lb=[]; uf_ub=[];
            for k=1:length(obj.controls)
                u_lb = [u_lb(:); obj.controls(k).lb(:)];
                u_ub = [u_ub(:); obj.controls(k).ub(:)];
                u0_lb = [u0_lb(:); obj.controls(k).lb0(:)];
                u0_ub = [u0_ub(:); obj.controls(k).ub0(:)];
                uf_lb = [uf_lb(:); obj.controls(k).lbf(:)];
                uf_ub = [uf_ub(:); obj.controls(k).ubf(:)];
            end
            obj.u0_lb = u0_lb;
            obj.u0_ub = u0_ub;
            obj.u_lb  = u_lb;
            obj.u_ub  = u_ub;
            obj.uf_lb = uf_lb;
            obj.uf_ub = uf_ub;
        end
        
        function obj = vectorize_parameters(obj)
            obj.p = obj.parameters.mx_vec();
            p_lb=[]; p_ub=[];
            for k=1:length(obj.parameters)
                p_lb = [p_lb(:); obj.parameters(k).lb(:)];
                p_ub = [p_ub(:); obj.parameters(k).ub(:)];
            end
            obj.p_lb  = p_lb;
            obj.p_ub  = p_ub;
        end
        
        function [t0, tf, t, x, u, p] = mx_vecs(obj)
            t0 = obj.independent_initial.mx_vec();
            tf = obj.independent_final.mx_vec();
            t = obj.independent.mx_vec();
            x = obj.states.mx_vec();
            u = obj.controls.mx_vec();
            p = obj.parameters.mx_vec();
        end
        
    end
    
    %% Transcription
    methods
        
        function n = nx(obj)
            n = numel(obj.x);
        end
        
        function n = nu(obj)
            n = numel(obj.u);
        end
        
        function n = np(obj)
            n = numel(obj.p);
        end
        
        function [t, x, u, p] = get_paramlist(obj)
            t = obj.t;
            x = obj.x;
            u = obj.u;
            p = obj.p;
        end
        
        function fn = expr_fn(obj, expr)
            obj.set_mx();
            e = fw_eval(expr);
            fn = casadi.Function('fn', {obj.t, obj.x, obj.u, obj.p}, {e});
        end
        
        function [t0_lb, t0_ub, tf_lb, tf_ub] = t_bd(obj)
            t0_lb=[]; t0_ub=[]; tf_lb=[]; tf_ub=[];
            for k=1:length(obj.independent_initial)
                t0_lb = [t0_lb(:); obj.independent_initial(k).lb(:)];
                t0_ub = [t0_ub(:); obj.independent_initial(k).ub(:)];
            end
            for k=1:length(obj.independent_final)
                tf_lb = [tf_lb(:); obj.independent_final(k).lb(:)];
                tf_ub = [tf_ub(:); obj.independent_final(k).ub(:)];
            end
        end
        
        function [p_lb, p_ub] = p_bd(obj)
            p_lb=[]; p_ub=[];
            for k=1:length(obj.parameters)
                p_lb = [p_lb(:); obj.parameters(k).lb(:)];
                p_ub = [p_ub(:); obj.parameters(k).ub(:)];
            end
        end
        
        function [bool, T] = fixed_horizon(obj)
            if obj.t0_lb==obj.t0_ub && obj.tf_lb==obj.tf_ub && all( ...
                    ~isinf( [obj.t0_lb; obj.t0_ub; obj.tf_lb; obj.tf_ub] ))
                bool = true;
            else
                bool = false;
            end    
            T = obj.tf_lb - obj.t0_lb;
        end
        
        function obj = set_mx(obj)
            for v=obj.variables()
                v.set_mx();
            end
        end 
        
        function obj = set_sym(obj)
            for v=obj.variables()
                v.set_sym();
            end
        end
        
        function obj = reset_sym(obj)
            for v=obj.variables()
                v.reset_value();
            end
        end
        
    end
        
    %% Presentation
    methods
        
        function present(obj)
            
            obj.independent.set_sym();
            obj.independent_initial.set_sym();
            obj.independent_final.set_sym();
            obj.states.set_sym();
            obj.algebraics.set_sym();
            obj.controls.set_sym();
            obj.parameters.set_sym();
            
            % objective function
            val = fw_eval(obj.objective);
            if isempty(obj.name)
                title = 'Optimal Control Problem';
            else
                title = obj.name;
            end
            fprintf(['[Yop] ', title, '\n']);
            fprintf('  min\t');
            if isnumeric(val)
                fprintf(num2str(val));
            else
                fprintf(char(val));
            end
            fprintf('\n');
            
            % Constraints
            fprintf('  s.t.\n');
            
            obj.print_box('independent_initial');
            obj.print_box('independent_final');
            
            if ~isempty(obj.states) || ~isempty(obj.controls)
                fprintf('  @t0\n');
            end
            obj.print_box_timed('states', 'lb0', 'ub0', '(t0)');
            obj.print_box_timed('controls', 'lb0', 'ub0', '(t0)');
            
            if ~isempty(obj.states) || ~isempty(obj.algebraics) || ~isempty(obj.controls)
                fprintf('  @t\n');
            end
            obj.print_box('states');
            obj.print_box('algebraics');
            obj.print_box('controls');
            
            if ~isempty(obj.states) || ~isempty(obj.controls)
                fprintf('  @tf\n');
            end
            obj.print_box_timed('states', 'lbf', 'ubf', '(tf)');
            obj.print_box_timed('controls', 'lbf', 'ubf', '(tf)');
            
            if ~isempty(obj.parameters)
                fprintf('  Parameters\n');
            end
            obj.print_box('parameters');
            
            
            if ~isempty(obj.differetial_eqs)
                fprintf('  ODE\n');
            end
            for k=1:length(obj.differetial_eqs)
                fprintf('\tder(');                                
                fprintf(char(fw_eval(obj.differetial_eqs(k).var)));
                fprintf(') == ');
                rhs = fw_eval(obj.differetial_eqs(k).expr);
                if length(char(rhs)) > 40
                    vs = yop.get_vars(obj.differetial_eqs(k).expr);
                    args = '';
                    for n=1:length(vs)
                        args = [args(:).', vs{n}.name, ', '];
                    end
                    fprintf(['f(', args(1:end-2), ')']);
                else
                    if isnumeric(rhs)
                        fprintf(num2str(rhs));
                    else
                        fprintf(char(rhs));
                    end
                end
                fprintf('\n');
            end
            
            if ~isempty(obj.eq)
                fprintf('  Equality\n');
            end
            for k=1:length(obj.eq)
                fprintf('\t');
                fprintf(char(fw_eval(obj.eq{k})));
                fprintf('\n');
            end
            
            if ~isempty(obj.ieq)
                fprintf('  Inequality\n');
            end
            for k=1:length(obj.ieq)
                fprintf('\t');
                fprintf(char(fw_eval(obj.ieq{k})));
                fprintf('\n');
            end
            
        end
        
        function obj = print_box(obj, varstr)
            % Should replace this implementation with box_timed and then
            % remove box_timed as function name.
            for v=obj.(varstr)
                var_name = v.var.name;
                ub = char(sym(v.ub));
                lb = char(sym(v.lb));
                bc = ['\t', lb, ' <= ', var_name, ' <= ', ub, '\n'];
                fprintf(bc);
            end
        end
        
        function print_box_timed(obj, varstr, lbstr, ubstr, timestr)
            for v=obj.(varstr)
                var_name = [v.var.name, timestr];
                ub = char(sym(v.(ubstr)));
                lb = char(sym(v.(lbstr)));
                bc = ['\t', lb, ' <= ', var_name, ' <= ', ub, '\n'];
                fprintf(bc);
            end
        end
        
    end
end