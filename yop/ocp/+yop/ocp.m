classdef ocp < handle
    properties
        % Misc
        name
        
        % objective
        objective
        
        % Variables
        independent         = yop.ocp_var.empty(1,0);
        independent_initial = yop.ocp_var.empty(1,0);
        independent_final   = yop.ocp_var.empty(1,0);
        states              = yop.ocp_var.empty(1,0);
        algebraics          = yop.ocp_var.empty(1,0);
        controls            = yop.ocp_var.empty(1,0);
        parameters          = yop.ocp_var.empty(1,0);
        differetial_eqs
        constraints
        timepoints
        integrals
        vars
        
        % Vectorized form
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
        
        ode
        
        eq
        ieq
        
        e2i
        i2e
    end
    
    properties (Hidden, Access=private)
        % Pertmutation of the state vector:
        %e2i        % Go from the order in which the states are found to how 
                   % the are orderered in the dynamics. This is the 
                   % internal order, becuase then it makes senso to do
                   % dx + x, which is used by the explicit rk integrators.
                   
        %i2e        % Go from how they are ordered in the derivative vector
                   % to the order in which the states are found. This is
                   % the order that is obtained by the formulation of the
                   % problem and is generally random. Here it is possible
                   % that the state variables are found in a different
                   % order than their derivatives, so it does not makes
                   % sense to use this order internally.
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
            obj.objective = yop.ocp_expr(objective);
        end
        
        function obj = max(obj, objective)
            % Maybe its better to note that it is a max problem, so that
            % the present function can present it as a maximization problem
            obj.objective = yop.ocp_expr(-objective);
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
            obj.find_variables_timepoints_integrals();
            obj.classify_variables();
            
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
            
            % 7) vectorize ocp
            obj.vectorize_variables();
            
            % 8) Compute timepoint and integral functions
            obj.set_timepoint_placeholders();
            obj.set_integral_placeholders();
            obj.set_timepoint_functions(); % måste komma efter vectorize 
            obj.set_integral_functions();  % för att parameterlistan ska vara rätt
            
            % 9) Convert to transcription invariant/variant form
            tri = yop.tr_invar(nbc); 
            obj.eq = tri.eq;
            obj.ieq = tri.ieq;
            
            % 10) Compute function objects for objective and  
            %     path constraints
            obj.compute_objective_function();
            obj.compute_pathcon_functions();
            
            % 11) Error checking
            
            % OCP built!
        end
        
        function obj = find_variables_timepoints_integrals(obj)
            obj.vars = {};
            obj.timepoints = yop.ocp_timepoint.empty(1,0);
            obj.integrals = yop.ocp_int.empty(1,0);
            
            exprs = {obj.objective.expr, obj.constraints{:}};
            
            % Hoist first iteration in order to avoid to visit the same 
            % node twice as all nodes are stored in 'visited', which is
            % reused.
            [tsort, n_elem, visited] = topological_sort(exprs{1});
            for n=1:n_elem
                obj.classify_node(tsort{n});
            end
            for k=2:length(exprs)
                [tsort, n_elem, visited] = ...
                    topological_sort(exprs{k}, visited);
                for n=1:n_elem
                    obj.classify_node(tsort{n});
                end
            end
        end
        
        function obj = classify_node(obj, node)
            if isa(node, 'yop.ast_variable')
                obj.vars = {obj.vars{:}, node};
            elseif isa(node, 'yop.ast_timepoint')
                obj.timepoints(end+1) = yop.ocp_timepoint(node);
            elseif isa(node, 'yop.ast_int')
                obj.integrals(end+1) = yop.ocp_int(node);
            end
        end
        
        function obj = classify_variables(obj)
            for k=1:length(obj.vars)
                vk = obj.vars{k};
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
            obj.states(end+1) = yop.ocp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            obj.algebraics(end+1) = yop.ocp_var(z);
        end
        
        function obj = add_control(obj, u)
            obj.controls(end+1) = yop.ocp_var(u);
        end
        
        function obj = add_parameter(obj, p)
            obj.parameters(end+1) = yop.ocp_var(p);
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
                        var.ub(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lb(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
                        var.ub0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lb0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
                        var.ubf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lbf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ubf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lbf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
                        var.ub(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lb(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
                        var.ub0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lb0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub0(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
                        var.ubf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                        var.lbf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lbf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ubf(re.idx_var) = ...
                            yop.get_subexpr(bnd, re.idx_expr);
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
        
        function obj = vectorize_variables(obj)
            obj.vectorize_independent();
            obj.vectorize_state();
            obj.vectorize_algebraic();
            obj.vectorize_control();
            obj.vectorize_parameters();
        end
        
        function obj = vectorize_independent(obj)
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
            %  that reason it is necessary to compute a perturbation
            %  vector, Px, so that the elements come in the same order as
            %  their derivatives.
            
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
            [re, nr, output_idx] = obj.reaching_states(ode_var);
            assert(all( diff(obj.states.get_enum()) == 1 ) ...
                && obj.states(1).enum(1)==1, '[Yop] Unexpected error.');
            
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
                    ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                    obj.differetial_eqs(end+1) = ...
                        yop.ocp_ode(vk, zeros(size(vk(:))));
                    output_idx = [output_idx(:); re(k).enum(idx)];
                end
            end
            % case iii) do not reach at all, so they are added in full.
            for k=1:length(nr)
                vk = nr(k).var;
                ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                obj.differetial_eqs(end+1) = ...
                        yop.ocp_ode(vk, zeros(size(vk(:))));
                output_idx = [output_idx(:); nr(k).enum(:)];
            end
            
            % 4) State perturbation vector
            [~, idx] = sort(output_idx);
            obj.e2i = output_idx;% ordering: Ord(der(x))==Ord(xx(P))
            obj.i2e = idx;       % ordering: Ord(der(x(P)))==Ord(xx)
            
            % -- Vectorize state --
            % 5) Set ode_rhs
            [tt, xx, uu, pp] = obj.mx_parameter_list();
            obj.set_mx();
            ode_fn = casadi.Function('x', {tt,xx,uu,pp}, ...
                {fw_eval(ode_expr)});
            obj.ode = casadi.Function('ode', {tt,xx,uu,pp}, ...
                {ode_fn(tt,xx(obj.i2e),uu,pp)});
            
            % 6) Set state bounds
            x0lb=[]; x0ub=[]; xlb=[]; xub=[]; xflb=[]; xfub=[];
            for k=1:length(obj.states)
                xlb = [xlb(:); obj.states(k).lb(:)];
                xub = [xub(:); obj.states(k).ub(:)];
                x0lb = [x0lb(:); obj.states(k).lb0(:)];
                x0ub = [x0ub(:); obj.states(k).ub0(:)];
                xflb = [xflb(:); obj.states(k).lbf(:)];
                xfub = [xfub(:); obj.states(k).ubf(:)];
            end
            obj.x0_lb = x0lb(obj.e2i);
            obj.x0_ub = x0ub(obj.e2i);
            obj.x_lb  = xlb(obj.e2i);
            obj.x_ub  = xub(obj.e2i);
            obj.xf_lb = xflb(obj.e2i);
            obj.xf_ub = xfub(obj.e2i);
        end        
        
        function obj = vectorize_algebraic(obj)
            zlb=[]; zub=[];
            for k=1:length(obj.algebraics)
                zlb = [zlb(:); obj.algebraics(k).lb(:)];
                zub = [zub(:); obj.algebraics(k).ub(:)];
            end
            obj.z_lb  = zlb;
            obj.z_ub  = zub;
        end
        
        function obj = vectorize_control(obj)
            u0lb=[]; u0ub=[]; ulb=[]; uub=[]; uflb=[]; ufub=[];
            for k=1:length(obj.controls)
                ulb = [ulb(:); obj.controls(k).lb(:)];
                uub = [uub(:); obj.controls(k).ub(:)];
                u0lb = [u0lb(:); obj.controls(k).lb0(:)];
                u0ub = [u0ub(:); obj.controls(k).ub0(:)];
                uflb = [uflb(:); obj.controls(k).lbf(:)];
                ufub = [ufub(:); obj.controls(k).ubf(:)];
            end
            obj.u0_lb = u0lb;
            obj.u0_ub = u0ub;
            obj.u_lb  = ulb;
            obj.u_ub  = uub;
            obj.uf_lb = uflb;
            obj.uf_ub = ufub;
        end
        
        function obj = vectorize_parameters(obj)
            plb=[]; pub=[];
            for k=1:length(obj.parameters)
                plb = [plb(:); obj.parameters(k).lb(:)];
                pub = [pub(:); obj.parameters(k).ub(:)];
            end
            obj.p_lb  = plb;
            obj.p_ub  = pub;
        end
        
        function obj = set_timepoint_placeholders(obj)
            k = 1;
           for tp=obj.timepoints
               if isa(tp.timepoint, 'yop.ast_independent_initial')
                   tp_text = 't0';
               elseif isa(tp.timepoint, 'yop.ast_independent_final')
                   tp_text = 'tf';
               else
                   tp_text = num2str(floor(tp.timepoint));% shouold relace . with _
               end
               nme = ['expr', num2str(k), '_t_eq_', tp_text]; 
               tp.mx= casadi.MX.sym(nme, size(tp.node,1), size(tp.node,2));
               tp.sym = sym(nme, size(tp.node));
               k = k+1;
           end
        end
        
        function obj = set_integral_placeholders(obj)
            k = 1;
           for tp=obj.integrals
               str = ['int', num2str(k)];
               tp.mx= casadi.MX.sym(str, size(tp.node,1), size(tp.node,2));
               tp.sym = sym(str, size(tp.node));
               k = k+1;
           end
        end
        
        function obj = set_timepoint_functions(obj)
            for tp = obj.timepoints
                tp.fn = obj.mx_function_object(tp.expr);
            end
        end
        
        function obj = set_integral_functions(obj)
            for int = obj.integrals
                int.fn = obj.mx_function_object(int.expr);
            end
        end
        
        function obj = compute_objective_function(obj)
            J = obj.objective; 
            J.fn = obj.mx_function_object(J.expr);
        end
        
        function obj = compute_pathcon_functions(obj)
            for pc = [obj.eq, obj.ieq]
                pc.fn = obj.mx_function_object(pc.expr);
            end
        end
        
    end
    
    %% Internal methods and helper methods
    
    methods (Access=private)
        
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
        
        function [ast_var, ocp_var] = find_variable(obj, id)
            for ocp_var = obj.variables
                if ocp_var.var.id == id
                    ast_var = ocp_var;
                    return;
                end
            end
            error('[Yop] Error: ID not found');
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
        
        function obj = set_mx(obj)
           obj.variables.set_mx();
           obj.timepoints.set_mx();
           obj.integrals.set_mx();
        end
        
        function [tt, xx, uu, pp] = mx_parameter_list(obj)
            % Used to derive function objects. Does not perturb state
            % vector, which is why this function is private. All outide
            % communication should have no knowledge of that the state
            % vector needs to be perurbed.
            tt = obj.independent.mx_vec();
            xx = obj.states.mx_vec();
            uu = obj.controls.mx_vec();
            pp = obj.parameters.mx_vec();
        end
        
        function fn = mx_function_object(obj, expr)
            obj.set_mx();
            [tt, xx, uu, pp] = obj.mx_parameter_list();
            tps = obj.timepoints.mx_vec();
            ints = obj.integrals.mx_vec();
            obj.set_mx();
            mx_expr = fw_eval(expr);
            innr = casadi.Function('i', {tt,xx,uu,pp,tps,ints}, {mx_expr});
            fn = casadi.Function('o', {tt,xx,uu,pp,tps,ints}, ...
                {innr(tt,xx(obj.i2e),uu,pp,tps,ints)});
        end
        
        function [re, nr, reaching_elements] = reaching_states(obj, expr)
            % Documentation for this function is that of yop.reaching_elems
            % as it is a variation on that one, with the difference that it
            % only does it for the 'states' property only. The reason for
            % this implementation is that computation of the perturbation
            % function for the state vector relies on that the states are
            % enumerated in the order they appear in the state vector.
            % Also, a sligth speed up is obtained by saving the enumerated
            % index to the states(k) property directly.
            xs = arrayfun(@(e) e.var, obj.states, 'UniformOutput', false);
            re(length(xs)) = yop.re_data();
            e0 = 1;
            for k=1:length(re)
                re(k).var = xs{k};
                [e0, idx] = re(k).enumerate(e0); 
                obj.states(k).enum = idx;
            end
            reaching_elements = propagate_value(expr);
            reaching_elements(~isa_variable(expr)) = -1;
            for k=1:length(re)
                re(k).set_expr_idx(reaching_elements);
            end
            reaching = arrayfun(@(v) ~isempty(v.reaching), re);            
            nr = re(~reaching);
            re = re( reaching);
        end
        
    end
    
    
    %% Transcription
    methods
        
        function n = n_x(obj)
            n = numel(obj.states.mx_vec);
        end
        
        function n = n_u(obj)
            n = numel(obj.controls.mx_vec);
        end
        
        function n = n_p(obj)
            n = numel(obj.parameters.mx_vec);
        end
        
        function n = n_tp(obj)
            n = numel(obj.timepoints.mx_vec);
        end
        
        function n = n_int(obj)
            n = numel(obj.integrals.mx_vec);
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
            val = propagate_value(obj.objective.expr);
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
            
            if ~isempty(obj.states) || ~isempty(obj.controls)
                fprintf('  @t\n');
            end
            obj.print_box('states');
            %obj.print_box('algebraics');
            obj.print_box('controls');
            
            if ~isempty(obj.states) || ~isempty(obj.controls)
                fprintf('  @tf\n');
            end
            obj.print_box_timed('states', 'lbf', 'ubf', '(tf)');
            obj.print_box_timed('controls', 'lbf', 'ubf', '(tf)');
            
            if ~(size(obj.parameters.var,1)==0 && ...
                    size(obj.parameters.var,2)==0)
                fprintf('  Parameters\n');
                obj.print_box('parameters');
            end
            
            
            if ~isempty(obj.differetial_eqs)
                fprintf('  ODE\n');
            end
            for k=1:length(obj.differetial_eqs)
                fprintf('\tder(');                                
                fprintf(char(propagate_value(obj.differetial_eqs(k).var)));
                fprintf(') == ');
                rhs = propagate_value(obj.differetial_eqs(k).expr);
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
%                 fprintf(char(propagate_value(obj.eq{k})));
                fprintf('\n');
            end
            
            if ~isempty(obj.ieq)
                fprintf('  Inequality\n');
            end
            for k=1:length(obj.ieq)
                fprintf('\t');
%                 fprintf(char(propagate_value(obj.ieq{k})));
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