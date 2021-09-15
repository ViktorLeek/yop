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
        odes
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
        % Permutation of the state vector:
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
            % Find all variables from the relations and expressions that
            % make up the problem
            obj.find_variables_timepoints_integrals();
            obj.classify_variables();
            
            % SRF
            srf = yop.ocp.to_srf(obj.constraints);
            [box, nbox] = yop.ocp.separate_box(srf);
            
            % Box
            [box_t, box_t0, box_tf] = yop.ocp.timed_box(box);
            obj.parse_box(box_t, 'lb', 'ub');
            obj.parse_box(box_t0, 'lb0', 'ub0');
            obj.parse_box(box_tf, 'lbf', 'ubf');
            obj.set_box_bounds(); % Includes processing default values
            
            % Sort the nonbox constraints: ode, inequality, equality
            [odes_, alg, eqs, ieqs] = yop.ocp.sort_nonbox(nbox);
            obj.odes = odes_;
            obj.eq = yop.ocp.split_transcription_invariance(eqs);
            obj.ieq = yop.ocp.split_transcription_invariance(ieqs);
            
            % Vectorize ocp
            obj.vectorize_ocp();
            
            % Error checking
            
            % OCP built!
        end
        
        function obj = find_variables_timepoints_integrals(obj)
            %______________________________________________________________
            %|YOP.OCP/FIND_VARIABLE_TIMEPOINTS_INTEGRALS From all problem |
            %|expressions finds all variables, timepoints and expressions.|
            %|                                                            |
            %| Use:                                                       |
            %|   find_variables_timepoints_integrals(obj)                 |
            %|   obj.find_variables_timepoints_integrals()                |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            obj.vars = {};
            obj.timepoints = yop.ocp_timepoint.empty(1,0);
            obj.integrals = yop.ocp_int.empty(1,0);
                        
            % Hoist first iteration in order to avoid to visit the same
            % node twice as all nodes are stored in 'visited', which is
            % reused.
            exprs = {obj.objective.expr, obj.constraints{:}};
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
            %______________________________________________________________
            %|YOP.OCP/CLASSIFY_NODE Classify the node depending on being  |
            %|a variable, timepoint, or integral.                         |
            %|                                                            |
            %| Use:                                                       |
            %|   classify_node(obj, node)                                 |
            %|   obj.classify_node(node)                                  |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   node - AST node to classify.                             |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            if isa(node, 'yop.ast_variable')
                obj.vars = {obj.vars{:}, node};
            elseif isa(node, 'yop.ast_timepoint')
                obj.timepoints(end+1) = yop.ocp_timepoint(node);
            elseif isa(node, 'yop.ast_int')
                obj.integrals(end+1) = yop.ocp_int(node);
            end
        end
        
        function obj = classify_variables(obj)
            %______________________________________________________________
            %|YOP.OCP/CLASSIFY_VARIABLE Classify variables based on       |
            %|independent, independent initial, independent final,        |
            %|state, algebraic, control input, free parameter.            |
            %|                                                            |
            %| Use:                                                       |
            %|   classify_variables(obj)                                  |
            %|   obj.classify_variables()                                 |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
            %______________________________________________________________
            %|YOP.OCP/ADD_INDEPENDENT Add independent variable to the OCP.|
            %|                                                            |
            %| Use:                                                       |
            %|   add_independent(obj, t)                                  |
            %|   obj.add_independent(t)                                   |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   t - The ast_independent.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            if isempty(obj.independent)
                obj.independent = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'independent variable']);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            %______________________________________________________________
            %|YOP.OCP/ADD_INDEPENDENT_INITIAL Add initial time parameter  |
            %|to the OCP.                                                 |
            %|                                                            |
            %| Use:                                                       |
            %|   add_independent_initial(obj, t)                          |
            %|   obj.add_independent_initial(t)                           |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   t - The ast_independent_initial.                         |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'initial bound for the independent variable']);
            end
        end
        
        function obj = add_independent_final(obj, t)
            %______________________________________________________________
            %|YOP.OCP/ADD_INDEPENDENT_FINAL Add end time parameter to the |
            %|OCP.                                                        |
            %|                                                            |
            %| Use:                                                       |
            %|   add_independent_final(obj, t)                            |
            %|   obj.add_independent_final(t)                             |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   t - The ast_independent_final.                           |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'terminal bound for the independent variable']);
            end
        end
        
        function obj = add_state(obj, x)
            %___________________________________________________
            %|YOP.OCP/ADD_STATE Add state variable to the OCP. |
            %|                                                 |
            %| Use:                                            |
            %|   add_state(obj, x)                             |
            %|   obj.add_state(x)                              |
            %|                                                 |
            %| Parameters:                                     |
            %|   obj - Handle to the ocp.                      |
            %|   x - The ast_state.                            |
            %|                                                 |
            %| Return values:                                  |
            %|   obj - Handle to the ocp.                      |
            %|_________________________________________________|
            obj.states(end+1) = yop.ocp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            %___________________________________________________________
            %|YOP.OCP/ADD_ALGEBRAIC Add algebraic variable to the OCP. |
            %|                                                         |
            %| Use:                                                    |
            %|   add_algebraic(obj, z)                                 |
            %|   obj.add_algebraic(z)                                  |
            %|                                                         |
            %| Parameters:                                             |
            %|   obj - Handle to the ocp.                              |
            %|   z - The ast_algebraic.                                |
            %|                                                         |
            %| Return values:                                          |
            %|   obj - Handle to the ocp.                              |
            %|_________________________________________________________|
            obj.algebraics(end+1) = yop.ocp_var(z);
        end
        
        function obj = add_control(obj, u)
            %____________________________________________________
            %|YOP.OCP/ADD_CONTROL Add control input to the OCP. |
            %|                                                  |
            %| Use:                                             |
            %|   add_control(obj, u)                            |
            %|   obj.add_control(u)                             |
            %|                                                  |
            %| Parameters:                                      |
            %|   obj - Handle to the ocp.                       |
            %|   u - The ast_control.                           |
            %|                                                  |
            %| Return values:                                   |
            %|   obj - Handle to the ocp.                       |
            %|__________________________________________________|
            obj.controls(end+1) = yop.ocp_var(u);
        end
        
        function obj = add_parameter(obj, p)
            %___________________________________________________________
            %|YOP.OCP/ADD_PARAMETER Add free optimization parameter to |
            %|the OCP.                                                 |
            %|                                                         |
            %| Use:                                                    |
            %|   add_parameter(obj, p)                                 |
            %|   obj.add_parameter(p)                                  |
            %|                                                         |
            %| Parameters:                                             |
            %|   obj - Handle to the ocp.                              |
            %|   p - The ast_parameter.                                |
            %|                                                         |
            %| Return values:                                          |
            %|   obj - Handle to the ocp.                              |
            %|_________________________________________________________|
            obj.parameters(end+1) = yop.ocp_var(p);
        end
        
        function obj = set_box_bounds(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_BOX_BOUNDS From box constraints and default     |
            %|values, set the box constraints for the problem.            |
            %|                                                            |
            %| Use:                                                       |
            %|   set_box_bounds(obj)                                      |
            %|   obj.set_box_bounds()                                     |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
        
        function obj = vectorize_ocp(obj)
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_OCP Vectorize the OCP                     |
            %|                                                            |
            %| Use:                                                       |
            %|   vectorize_ode(obj)                                       |
            %|   obj.vectorize_ode()                                      |
            %|                                                            |
            %| Description:                                               |
            %|   Vectorizes the OCP. The imporant thing to keep in mind   |
            %|   is that the ode is vectorized before any mx function     |
            %|   expression is created. The reason is that it computes    |
            %|   the state permutation vectors, which are necessar for    |
            %|   providing a good interface for other functions and       |
            %|   classes. Because of this, vectorize_ode is run first.    |
            %|   Also, timepoint and integral placeholders must be set    |
            %|   before their functions can be created.                   |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            obj.vectorize_ode();
            obj.vectorize_independent();
            obj.vectorize_state();
            obj.vectorize_algebraic();
            obj.vectorize_control();
            obj.vectorize_parameters();
            obj.set_timepoint_placeholders();
            obj.set_integral_placeholders();
            obj.set_timepoint_functions();
            obj.set_integral_functions();
            obj.compute_objective_function();
            obj.compute_pathcon_functions();
        end
        
        function obj = vectorize_ode(obj)
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_ODE Vectorize all ODEs and set default    |
            %|odes for states that are not bound by an ode.               |
            %|                                                            |
            %| Use:                                                       |
            %|   vectorize_ode(obj)                                       |
            %|   obj.vectorize_ode()                                      |
            %|                                                            |
            %| Description:                                               |
            %|   Analyzes which of the state elements the reaches the ode |
            %|   rhs: der(var) == expr. The elements of the entire state  |
            %|   vector that does not reach is given the default          |
            %|   derivative: der(var) == 0.                               |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
  
            
            % ode_var - Variables of ode rhs (ode(var) == expr) in
            %           vectorized form.
            % ode_expr - Expressions of ode lhs (ode(var) == expr) in
            %            vectorized form.
            % The ode is canonicalized, so all variables are to the left
            % and expressions to the right.
            ode_var = []; ode_expr = [];
            for k=1:length(obj.odes)
                ode_var = [ode_var(:); obj.odes{k}.lhs(:)];
                ode_expr = [ode_expr(:); obj.odes{k}.rhs(:)];
            end
            
            % Analyze the elements that reaches the ode rhs: ode(var)==...
            % The purpose is to set the ode of those that does not reach to
            % zero. re - reaching elements, nr - not reaching
            [re, nr] = obj.reaching_states(ode_var);
            
            % Test if some of the elements do not reach
            for k=1:length(re)
                if numel(re(k).enum) ~= numel(re(k).reaching)
                    [~, idx] = setdiff(re(k).enum, re(k).reaching);
                    vk = re(k).var(idx);
                    ode_var = [ode_var(:); vk(:)];
                    ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                    
                    if yop.settings.warnings
                        warning(yop.msg.default_der(vk.name));
                    end
                end
            end
            
            % Entire variable do not reach, set default ode for all
            for k=1:length(nr)
                vk = nr(k).var;
                ode_var = [ode_var(:); vk(:)];
                ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                
                if yop.settings.warnings
                    warning(yop.msg.default_der(vk.name));
                end
            end
            
            % Very important step! Necessary before any MX function object
            % can be set.
            obj.compute_permutation_vector(ode_var);
            
            % Compute ODE function in two steps
            obj.ode = obj.mx_ode_function(ode_expr);
        end
        
        function obj = compute_permutation_vector(obj, ode_var)
            %______________________________________________________________
            %|YOP.OCP/COMPUTE_PERMUTATION_VECTOR Computes the vector that |
            %|permutes the state vector in such a way that the vector     |
            %|elemnts appear in the same order as they to in the state    |
            %|derivative.                                                 |
            %|                                                            |
            %| Use:                                                       |
            %|   compute_permutation_vector(obj, ode_var)                 |
            %|   obj.compute_permutation_vector(ode_var)                  |
            %|                                                            |
            %| Description:                                               |
            %|   Because the integration methods might add state and      |
            %|   deriative in order to compute the next step, it is       |
            %|   important that elements appear in the same order in the  |
            %|   state variable vector as in the state derivative vector. |
            %|   Yop solves that by computing two permutation vectors.    |
            %|   One maps from external to internal representation (e2i)  |
            %|   and one maps from internal to external representation    |
            %|   (i2e). External representation could be any order and is |
            %|   simply the order in which the states are detected.       |
            %|   Internal order is also random to some degree since it    |
            %|   might depend on the order in which odes are detected.    |
            %|   Nevertheless, it is necessary to work with only one      |
            %|   representation internally, which is why yop computes     |
            %|   these permutation vectors.                               |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   ode_var - AST expression for the variables as they       |
            %|             appear in the ODE.                             |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            [~, ~, output_idx] = obj.reaching_states(ode_var);
            [~, idx] = sort(output_idx);
            obj.e2i = output_idx;
            obj.i2e = idx;
        end
        
        function obj = vectorize_independent(obj)
            obj.t0_lb = obj.independent_initial.lb;
            obj.t0_ub = obj.independent_initial.ub;
            obj.tf_lb = obj.independent_final.lb;
            obj.tf_ub = obj.independent_final.ub;
        end
        
        function obj = vectorize_state(obj)
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_STATE Vectorize the OCP states.           |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_ALGEBRAIC Vectorize the OCP algebraic     |
            %|variables.                                                  |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            zlb=[]; zub=[];
            for k=1:length(obj.algebraics)
                zlb = [zlb(:); obj.algebraics(k).lb(:)];
                zub = [zub(:); obj.algebraics(k).ub(:)];
            end
            obj.z_lb  = zlb;
            obj.z_ub  = zub;
        end
        
        function obj = vectorize_control(obj)
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_CONTROL Vectorize the OCP control inputs. |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_PARAMETERS Vectorize the OCP free         |
            %|parameters.                                                 |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            plb=[]; pub=[];
            for k=1:length(obj.parameters)
                plb = [plb(:); obj.parameters(k).lb(:)];
                pub = [pub(:); obj.parameters(k).ub(:)];
            end
            obj.p_lb  = plb;
            obj.p_ub  = pub;
        end
        
        function obj = set_timepoint_placeholders(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_TIMEPOINT_PLACEHOLDERS Give every timepoint an  |
            %|MX variable of the same size as the expression to represent |
            %|its value.                                                  |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            k = 1;
            for tp=obj.timepoints
                if isa(tp.timepoint, 'yop.ast_independent_initial')
                    tp_text = 't0';
                elseif isa(tp.timepoint, 'yop.ast_independent_final')
                    tp_text = 'tf';
                else
                    tp_text = num2str(floor(tp.timepoint));
                end
                nme = ['expr', num2str(k), '_t_eq_', tp_text];
                tp.mx= casadi.MX.sym(nme, size(tp.node,1), size(tp.node,2));
                tp.sym = sym(nme, size(tp.node));
                k = k+1;
            end
        end
        
        function obj = set_integral_placeholders(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_INTEGRAL_PLACEHOLDERS Give every integral an MX |
            %|variable of the same size as the expression to represent    |
            %|its value.                                                  |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            k = 1;
            for tp=obj.integrals
                str = ['int', num2str(k)];
                tp.mx= casadi.MX.sym(str, size(tp.node,1), size(tp.node,2));
                tp.sym = sym(str, size(tp.node));
                k = k+1;
            end
        end
        
        function obj = set_timepoint_functions(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_TIMEPOINT_FUNCTIONS Compute mx function         |
            %|objects for the timepoint expressions.                      |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            for tp = obj.timepoints
                tp.fn = obj.mx_function_object(tp.expr);
            end
        end
        
        function obj = set_integral_functions(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_INTEGRAL_FUNCTIONS Compute mx function          |
            %|objects for the integral expressions.                       |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            for int = obj.integrals
                int.fn = obj.mx_function_object(int.expr);
            end
        end
        
        function obj = compute_objective_function(obj)
            %______________________________________________________________
            %|YOP.OCP/COMPUTE_OBJECTIVE_FUNCTION Compute mx function      |
            %|object for the objective function.                          |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            J = obj.objective;
            J.fn = obj.mx_function_object(J.expr);
        end
        
        function obj = compute_pathcon_functions(obj)
            %______________________________________________________________
            %|YOP.OCP/COMPUTE_PATHCON_FUNCTION Compute mx function        |
            %|objects for the path constraints.                           |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            for pc = [obj.eq, obj.ieq]
                pc.fn = obj.mx_function_object(pc.expr);
            end
        end
        
    end
    
    %% Internal methods and helper methods
    
    methods (Access=private)
        
        function obj = set_box_bnd(obj, var_field, lb_field, ...
                ub_field, lb_def, ub_def)
            %______________________________________________________________
            %|YOP.OCP/SET_BOX_BND Set bound on a problem variable.        |
            %|                                                            |
            %| Use:                                                       |
            %|  set_box_bnd(obj,var_field,lb_field,ub_field,lb_def,ub_def)|
            %|  obj.set_box_bnd(var_field,lb_field,ub_field,lb_def,ub_def)|
            %|                                                            |
            %| Example:                                                   |
            %|   obj.set_box_bnd('states', 'lb', 'ub', -inf, inf)         |
            %|   obj.set_box_bnd('states', 'lb0', 'ub0', 0, 0)            |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   var_field - String for the variable field/property.      |
            %|   lb_field - String for the lower bound field/property.    |
            %|   ub_field - String for the upper bound field/property.    |
            %|   lb_def - The default lower bound value.                  |
            %|   ub_def - The default upper bound value.                  |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
            %______________________________________________________________
            %|YOP.OCP/FIND_VARIABLE Find a problem variable by its ID.    |
            %|Error if the ID is not found among the problem variables.   |
            %|                                                            |
            %| Use:                                                       |
            %|   [ast_var, ocp_var] = find_variable(obj, id)              |
            %|   [ast_var, ocp_var] = obj.find_variable(id)               |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   id - Id of the sought variable.                          |
            %|                                                            |
            %| Return values:                                             |
            %|   ast_var - The ast node of the variable.                  |
            %|   ocp_var - The ocp information on the variable.           |
            %|____________________________________________________________|
            for ocp_var = obj.variables
                if ocp_var.var.id == id
                    ast_var = ocp_var;
                    return;
                end
            end
            error('[Yop] Error: ID not found');
        end
        
        function vars = variables(obj)
            %______________________________________________________________
            %|YOP.OCP/VARIABLES Get all OCP problem variables as a vector.|
            %|                                                            |
            %| Use:                                                       |
            %|   vars = variables(obj)                                    |
            %|   vars = obj.variables()                                   |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   vars - All problem variables vectorized.                 |
            %|____________________________________________________________|
            vars = [obj.independent(:).', ...
                obj.independent_initial(:).', ...
                obj.independent_final(:).', ...
                obj.states(:).', ...
                obj.algebraics(:).', ...
                obj.controls(:).', ...
                obj.parameters(:).'];
        end
        
        function obj = set_mx(obj)
            %______________________________________________________________
            %|YOP.OCP/SET_MX Set all OCP variables to their mx value.     |
            %|                                                            |
            %| Use:                                                       |
            %|   set_mx(obj)                                              |
            %|   obj.set_mx()                                             |
            %|                                                            |
            %| Description:                                               |
            %|   In order to get MX expressions for AST expressions it is |
            %|   necessary to set the value of all expression variables   |
            %|   to MX, which is done by this function.                   |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            obj.variables.set_mx();
            obj.timepoints.set_mx();
            obj.integrals.set_mx();
        end
        
        function [tt, xx, uu, pp] = mx_parameter_list(obj)
            %______________________________________________________________
            %|YOP.OCP/MX_PARAMETER_LIST Get external MX form of the       |
            %|expression parameter list. The function is private since    |
            %|the state vector is in external ordering which easily leads |
            %|to errors if not cautious.                                  |
            %|                                                            |
            %| Use:                                                       |
            %|   [tt, xx, uu, pp] = mx_parameter_list(obj)                |
            %|   [tt, xx, uu, pp] = obj.mx_parameter_list()               |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|                                                            |
            %| Return values:                                             |
            %|   tt - MX scalar for the independent variable.             |
            %|   xx - MX vector for the state variable.                   |
            %|   uu - MX vector for the control input.                    |
            %|   pp - MX vector for the free parameters.                  |
            %|____________________________________________________________|
            tt = obj.independent.mx_vec();
            xx = obj.states.mx_vec();
            uu = obj.controls.mx_vec();
            pp = obj.parameters.mx_vec();
        end
        
        function fn = mx_function_object(obj, expr)
            %______________________________________________________________
            %|YOP.OCP/MX_FUNCTION_OBJECT Compute function object from     |
            %|expression with the parameter list (t,x,u,p,tps,ints).      |
            %|                                                            |
            %| Use:                                                       |
            %|   fn = mx_ode_function(obj, expr)                          |
            %|   fn = obj.mx_ode_function(expr)                           |
            %|                                                            |
            %| Description:                                               |
            %|   Computes a casadi function object from an ast expression.|
            %|   The function object has its state parameter on the       |
            %|   internal form which means that                           |
            %|                                                            |
            %|       ordering(der(x)) == ordering(expr)                   |
            %|                                                            |
            %|   as opposed to external form where the states might be    |
            %|   permuted so that                                         |
            %|                                                            |
            %|       ordering(der(x)) != ordering(expr)                   |
            %|                                                            |
            %|   Function object parameters are:                          |
            %|     t - Independent variable.                              |
            %|     x - State vector (internal ordering).                  |
            %|     u - Control input.                                     |
            %|     p - Free parameters.                                   |
            %|     tps - Timepoint expressions.                           |
            %|     ints - Integral expressions.                           |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   expr - AST of expression to create a function object     |
            %|          from.                                             |
            %|                                                            |
            %| Return values:                                             |
            %|   fn - A function object when called with the arguments    |
            %|        t, x, u, p, tps, int (in that exact order, empty    |
            %|        ones are set to empty: []) computes 'expr'.         |
            %|____________________________________________________________|
            [tt, xx, uu, pp] = obj.mx_parameter_list();
            tps = obj.timepoints.mx_vec();
            ints = obj.integrals.mx_vec();
            args = {tt,xx,uu,pp,tps,ints};
            
            % Since we are seeking an MX expression, all ast_variable in
            % the ocp set their values to MX. Otherwise fw_eval could
            % result in anything.
            obj.set_mx();
            
            % First a function object is created that maps the expression
            % to variables in the order they appear in the state variable
            % vector. The state element order here is external, as states
            % simply appear as they were found. An idea is of course to
            % change the order of the elements in the variable vector, but
            % that is a bit cumbersome as one needs to keep track of the
            % individual elements, not just the state variables themselves.
            inner_expr = {fw_eval(expr)};
            innr_fn = casadi.Function('i', args, inner_expr);
            
            % To provide the internal interface the argument order is
            % changed. Here 'args' is reused. Here it simply maps inputs to
            % parameters in the function call. It could just as well have
            % been other mx variables with the same dimensions, so do not
            % be confused by the reuse of variable.
            outer_expr = {innr_fn(tt, xx(obj.i2e), uu, pp, tps, ints)};
            fn = casadi.Function('o', args, outer_expr);
        end
        
        function fn = mx_ode_function(obj, expr)
            %______________________________________________________________
            %|YOP.OCP/MX_ODE_FUNCTION Compute function object from        |
            %|expression with the ode parameter list (t,x,u,p).           |
            %|                                                            |
            %| Use:                                                       |
            %|   fn = mx_ode_function(obj, expr)                          |
            %|   fn = obj.mx_ode_function(expr)                           |
            %|                                                            |
            %| Description:                                               |
            %|   Computes a casadi function object from an ast expression.|
            %|   The function object has its state parameter on the       |
            %|   internal form which means that                           |
            %|                                                            |
            %|       ordering(der(x)) == ordering(expr)                   |
            %|                                                            |
            %|   as opposed to external form where the states might be    |
            %|   permuted so that                                         |
            %|                                                            |
            %|       ordering(der(x)) != ordering(expr)                   |
            %|                                                            |
            %|   Function object parameters are:                          |
            %|     t - Independent variable.                              |
            %|     x - State vector (internal ordering).                  |
            %|     u - Control input.                                     |
            %|     p - Free parameters.                                   |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   expr - AST of expression to create a function object     |
            %|          from.                                             |
            %|                                                            |
            %| Return values:                                             |
            %|   fn - A function object when called with the arguments    |
            %|        t, x, u, p (in that exact order, empty ones are set |
            %|        to empty: []) computes 'expr'.                      |
            %|____________________________________________________________|
            [tt, xx, uu, pp] = obj.mx_parameter_list();
            args = {tt,xx,uu,pp};
            
            % Since we are seeking an MX expression, all ast_variable in
            % the ocp set their values to MX. Otherwise fw_eval could
            % result in anything.
            obj.set_mx();
            
            % First a function object is created that maps the expression
            % to variables in the order they appear in the state variable
            % vector. The state element order here is external, as states
            % simply appear as they were found. An idea is of course to
            % change the order of the elements in the variable vector, but
            % that is a bit cumbersome as one needs to keep track of the
            % individual elements, not just the state variables themselves.
            inner_expr = {fw_eval(expr)};
            innr_fn = casadi.Function('i', args, inner_expr);
            
            % To provide the internal interface the argument order is
            % changed. Here 'args' is reused. Here it simply maps inputs to
            % parameters in the function call. It could just as well have
            % been other mx variables with the same dimensions, so do not
            % be confused by the reuse of variable.
            outer_expr = {innr_fn(tt, xx(obj.i2e), uu, pp)};
            fn = casadi.Function('o', args, outer_expr);
        end
        
        function [re, nr, reaching_enumeration] = reaching_states(obj, expr)
            %______________________________________________________________
            %|YOP.OCP/REACHING_STATES Reaching elements analysis for the  |
            %|OCP states given an expression on which to do the analysis. |
            %|                                                            |
            %| Use:                                                       |
            %|   [re,nr,reaching_enumeration] = reaching_states(obj, expr)|
            %|   [re,nr,reaching_enumeration] = obj.reaching_states(expr) |
            %|                                                            |
            %| Description:                                               |
            %|   Determines which of all state elements reaches expr.     |
            %|   Does the analysis by enumerating all state elements and  |
            %|   propagating them through the expression. Those elements  |
            %|   that are not state variables are set to zero. The        |
            %|   benefit of this is that it is possible to add two        |
            %|   reaching elements analyses and the reaching element's    |
            %|   enumeration is unaffected by those that are not states.  |
            %|                                                            |
            %| Parameters:                                                |
            %|   expr - Ast node for the expression of interest.          |
            %|                                                            |
            %| Return values:                                             |
            %|   re - A vector of yop.re_data, one element for every      |
            %|        variable that is reaching the expression.           |
            %|   nr - A vector of yop.re_data, one element for every      |
            %|        variable that is not reaching the expression.       |
            %|   reaching_enumeration - The enumeration that is reaching  |
            %|                          the expression. Zeroes for those  |
            %|                          that are not variables.           |
            %|____________________________________________________________|
            
            % Cell array with all states (ocp_variable objects)
            xs = arrayfun(@(e) e.var, obj.states, 'UniformOutput', false);
            re(length(xs)) = yop.re_data();
            
            % Enumerate all variable elements
            e0 = 1; % Enumeration start value.
            for k=1:length(re)
                re(k).var = xs{k};
                [e0, idx] = re(k).enumerate(e0);
                obj.states(k).enum = idx;
            end
            
            % Propagate elements through the expression
            reaching_enumeration = propagate_value(expr);
            
            % Filter
            reaching_enumeration(~(isa_der(expr) & isa_state(expr))) = 0;
            
            % Store results
            for k=1:length(re)
                re(k).set_expr_elements_reached(reaching_enumeration);
            end
            
            % Split reaching and not reaching
            reaching = arrayfun(@(v) ~isempty(v.reaching), re);
            nr = re(~reaching);
            re = re( reaching);
        end
        
    end
    
    %% Passes
    methods
        
        function obj = parse_box(obj, box, lb, ub)
            %______________________________________________________________
            %|YOP.OCP/PARSE_BOX Parses the box constraints and saves the  |
            %|the information to the probem variables.                    |
            %|                                                            |
            %| Use:                                                       |
            %|   parse_box(obj, box, lb, ub)                              |
            %|   obj.parse_box(box, lb, ub)                               |
            %|                                                            |
            %| example:                                                   |
            %|   obj.parse_box(box, 'lb', 'ub')                           |
            %|   obj.parse_box(box, 'lb0', 'ub0')                         |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - A handle to the ocp.                               |
            %|   box - Cell array with box constraints. They must be      |
            %|         canonicalized and refer to a single variable.      |
            %|   lb - A string for the field to which the lower bound     |
            %|        apply. See yop.ocp_var for the fields.              |
            %|   ub - A string for the field to which the upper bound     |
            %|        apply. See yop.ocp_var for the fields.              |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - A handle to the ocp.                               |
            %|____________________________________________________________|
            
            for k=1:length(box)
                re = yop.reaching_elems(box{k}.lhs);
                bnd = yop.prop_num(box{k}.rhs);
                var = obj.find_variable(re.var.id);
                switch class(box{k})
                    case 'yop.ast_eq' % var == bnd
                        var.(ub)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                        var.(lb)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.(ub)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.(lb)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
        end
        
    end
    
    methods (Static)
        
        function srf = to_srf(constraints)
            %______________________________________________________________
            %|YOP.OCP.TO_SRF To single relation form.                     |
            %|                                                            |
            %| Use:                                                       |
            %|   yop.ocp.to_srf(constraints)                              |
            %|   srf = yop.ocp.to_srf(constraints)                        |
            %|                                                            |
            %| Description:                                               |
            %|   A pass for converting all constraints to a form where    |
            %|   there is only one relation per constraint. For instance, |
            %|   the constraint expr1 <= expr2 <= expr3 is converted into |
            %|      i) expr1 <= expr2                                     |
            %|     ii) expr2 <= expr3                                     |
            %|   The purspose of this pass is to be simplify the job for  |
            %|   the remaning passes as they can start from the assumption|
            %|   that there is only one relation per constraint.          |
            %|                                                            |
            %| Parameters:                                                |
            %|   constraints - A cell array with constraints.             |
            %|                                                            |
            %| Return values:                                             |
            %|   srf - A cell array with the constraints on srf.          |
            %|____________________________________________________________|
            srf = {};
            for n=1:length(constraints)
                cn = constraints{n};
                relations = get_relations(cn);
                for k=1:length(relations)
                    rk = relations{k};
                    fnh = get_constructor(rk);
                    srf{end+1} = fnh(rmost(rk.lhs), lmost(rk.rhs));
                end
            end
        end
        
        function [box, nbox] = separate_box(srf)
            %______________________________________________________________
            %|YOP.OCP.SEPARATE_BOX Separates box from non-box constraints.|
            %|                                                            |
            %| Use:                                                       |
            %|   [box, nbox] = yop.ocp.separate_box(srf)                  |
            %|                                                            |
            %| Description:                                               |
            %|   Searches the constraints on srf form for constraints     |
            %|   that are box constraints. The requirements are dictated  |
            %|   by the yop.ocp.isa_box function.                         |
            %|                                                            |
            %| Parameters:                                                |
            %|   box - Cell array with canonicalized box constraints.     |
            %|   nbox - Cell array with non-box constraints.              |
            %|                                                            |
            %| Return values:                                             |
            %|   srf - Cell array containing constraints on srf form      |
            %|____________________________________________________________|
            box = {};
            nbox = {};
            for n=1:length(srf)
                var_num = yop.ocp.isa_box(srf{n}.lhs, srf{n}.rhs);
                num_var = yop.ocp.isa_box(srf{n}.rhs, srf{n}.lhs);
                isbox = var_num | num_var;
                box{end+1} = yop.get_subrel(srf{n}, isbox);
                nbox{end+1} = yop.get_subrel(srf{n}, ~isbox);
            end
            nbox = nbox(~cellfun('isempty', nbox));
            
            % Canonicalize
            box = box(~cellfun('isempty', box));
            for k=1:length(box)
                % Canonicalization is handled in the ast_relation node
                % since it is handled differently depending on class.
                box{k} = canonicalize_box(box{k});
            end
            
            % Split so that variables do not mix
            ubox = yop.ocp.unique_box(box);
            box = ubox;
        end
        
        function boolv = isa_box(var_cand, num_cand)
            %______________________________________________________________
            %|YOP.OCP.ISA_BOX Tests if the variable candidate and numeric |
            %|candidate form a box constraint.                            |
            %|                                                            |
            %| Use:                                                       |
            %|   boolv = yop.ocp.isa_box(var_cand, num_can)               |
            %|                                                            |
            %| Description:                                               |
            %|   Box constraints are made up of a variable and a numeric  |
            %|   bound. The function determines if var_cand is a variable |
            %|   and num_cand is a numeric value.                         |
            %|                                                            |
            %| Parameters:                                                |
            %|   var_cand - The expression for the variable candidate     |
            %|   num_cand - The expression for the numeric candidate      |
            %|                                                            |
            %| Return values:                                             |
            %|   boolv - A bool vector. If only one argument is a scalar, |
            %|           the bool vector has the same dimensions as athe  |
            %|           vector.                                          |
            %|____________________________________________________________|
            t0 = yop.initial_timepoint();
            tf = yop.final_timepoint();
            isvar = isa_variable(var_cand);
            isder = isa_der(var_cand);
            [istp, tps] = isa_timepoint(var_cand);
            isnum = isa_numeric(num_cand);
            boolv = isvar & ~isder & (~istp|tps==t0|tps==tf) & isnum;
        end
        
        function ubox = unique_box(cbox)
            %______________________________________________________________
            %|YOP.OCP.UNIQUE_BOX Sorts the box constraints so that every  |
            %|constraint only refers to one variable. Requires            |
            %|canonicalized box constraints.                              |
            %|                                                            |
            %| Use:                                                       |
            %|   uid_box = yop.ocp.unique_box(cbox)                       |
            %|                                                            |
            %| Description:                                               |
            %|   Uses the unique ID every node has to sort the            |
            %|   constraints so that every box constraint only refers to  |
            %|   one variable.                                            |
            %|                                                            |
            %| Parameters:                                                |
            %|   cbox - Cell array with canonicalized box constraints.    |
            %|                                                            |
            %| Return values:                                             |
            %|   ubox - Cell array with canonicalized box constraints     |
            %|          with only one variable per constraint.            |
            %|____________________________________________________________|
            ubox = {};
            for k=1:length(cbox)
                [~, ID] = isa_variable(cbox{k}.lhs);
                for uID = yop.row_vec(unique(ID))
                    ubox{end+1} = yop.get_subrel(cbox{k}, ID==uID);
                end
            end
            ubox = ubox(~cellfun('isempty', ubox));
        end
        
        function [bc_t, bc_t0, bc_tf] = timed_box(box)
            %______________________________________________________________
            %|YOP.OCP.TIMED_BOX Sorts the box constraints based on the    |
            %|timepoint information they convey. There are three          |
            %|categories: entire horizon, initial timepoint, terminal     |
            %|timepoint.                                                  |
            %|                                                            |
            %| Use:                                                       |
            %|   [bc_t, bc_t0, bc_tf] = yop.ocp.timed_box(box)            |
            %|                                                            |
            %| Parameters:                                                |
            %|   box - Cell array with canonicalized box constraints.     |
            %|                                                            |
            %| Return values:                                             |
            %|   bc - Cell array with canonicalized box constraints       |
            %|        that apply to the entire problem horizon.           |
            %|   bc_t0 - Cell array with canonicalized box constraints    |
            %|           that apply to the initial timepoint.             |
            %|   bc_tf - Cell array with canonicalized box constraints    |
            %|           that apply to the final timepoint.               |
            %|____________________________________________________________|
            t0 = yop.initial_timepoint();
            tf = yop.final_timepoint();
            bc_t={}; bc_t0={}; bc_tf={};
            for k=1:length(box)
                [istp, tp] = isa_timepoint(box{k}.lhs);
                bc_t{end+1} = yop.get_subrel(box{k}, ~istp);
                bc_t0{end+1} = yop.get_subrel(box{k}, tp==t0);
                bc_tf{end+1} = yop.get_subrel(box{k}, tp==tf);
            end
            bc_t = bc_t(~cellfun('isempty', bc_t));
            bc_t0 = bc_t0(~cellfun('isempty', bc_t0));
            bc_tf = bc_tf(~cellfun('isempty', bc_tf));
        end
        
        function [ode, alg, eq, ieq] = sort_nonbox(nbox)
            %______________________________________________________________
            %|YOP.OCP.SORT_NONBOC Sorts the constraints that are not box  |
            %|constraints. It returns constraints in the categories       |
            %|ordinary differential equation (ODE), algebraic equation,   |
            %|equality constraints, and inequality constraints.           |
            %|                                                            |
            %| Use:                                                       |
            %|   [ode, alg, eq, ieq] = yop.ocp.sort_nonbox(nbox)          |
            %|                                                            |
            %| Parameters:                                                |
            %|   nbox - Cell array with non-box constraints.              |
            %|                                                            |
            %| Return values:                                             |
            %|   ode - Cell array with ordinary differential equations.   |
            %|   alg - Cell array with algebraic equations.               |
            %|   eq - Cell array with equality constraints.               |
            %|   ieq - Cell array with inequality constraints.            |
            %|____________________________________________________________|
            ode={};alg={};eq={};ieq={};
            for k=1:length(nbox)
                if is_alg(nbox{k})
                    alg{end+1} = nbox{k};
                    
                elseif isa(nbox{k}, 'yop.ast_eq')
                    [~, ode_k, eq_k] = isa_ode(nbox{k});
                    ode{end+1} = ode_k;
                    eq{end+1} = eq_k;
                    
                else
                    ieq{end+1} = canonicalize(nbox{k});
                    
                end
            end
            ode = ode(~cellfun('isempty', ode));
            alg = alg(~cellfun('isempty', alg));
            ieq = ieq(~cellfun('isempty', ieq));
            eq  =  eq(~cellfun('isempty', eq));
        end
        
        function trcon = split_transcription_invariance(con)
            %______________________________________________________________
            %|YOP.OCP.SPLIT_TRANSCRIPTION_INVARIANCE Split constraints    |
            %|that are both transcription variant and invariant into two  |
            %|separate constraints.                                       |
            %|                                                            |
            %| Use:                                                       |
            %|   trcon = yop.ocp.split_transcription_invariance(con)      |
            %|                                                            |
            %| Description:                                               |
            %|   A constraint is transcription invariant if its           |
            %|   transcription does not change with the independent       |
            %|   variable. The distinction is important because           |
            %|   transcription variant constraints are assigned to every  |
            %|   point of the discretization grid, as opposed to an       |
            %|   invariant constraint that only needs to be applied       |
            %|   once. To make this easier to the transcription method    |
            %|   this function makes it so that transcription method      |
            %|   only needs to make a binary decision as to wheteher it   |
            %|   invariant or not, never something in between.            |
            %|                                                            |
            %| Parameters:                                                |
            %|   con - Cell array with canonlicalized constraints non-box |
            %|         constraints.                                       |
            %|                                                            |
            %| Return values:                                             |
            %|   trcon - Cell array with the constraints separated based  |
            %|           on their transcription invariane.                |
            %|____________________________________________________________|
            trcon = {};
            for k=1:length(con)
                invariant = is_transcription_invariant(con{k}.lhs);
                if all(invariant) || all(~invariant)
                    trcon{end+1} = con{k};
                else
                    trcon{end+1} = yop.get_subrelation(con{k}, invariant);
                    trcon{end+1} = yop.get_subrelation(con{k}, ~invariant);
                end
            end
        end
        
    end
    
    
    %% Transcription
    methods
        
        function n = n_x(obj)
            %___________________________________________________
            %|YOP.OCP/N_X Get the number of states in the OCP. |
            %|                                                 |
            %| Use:                                            |
            %|   n = n_x(obj)                                  |
            %|   n = obj.n_x()                                 |
            %|   n = obj.n_x                                   |
            %|                                                 |
            %| Parameters:                                     |
            %|   obj - A handle to the ocp.                    |
            %|                                                 |
            %| Return values:                                  |
            %|   n - The number of states.                     |
            %|_________________________________________________|
            n = numel(obj.states.mx_vec);
        end
        
        function n = n_u(obj)
            %___________________________________________________________
            %|YOP.OCP/N_U Get the number of control inputs in the OCP. |
            %|                                                         |
            %| Use:                                                    |
            %|   n = n_u(obj)                                          |
            %|   n = obj.n_u()                                         |
            %|   n = obj.n_u                                           |
            %|                                                         |
            %| Parameters:                                             |
            %|   obj - A handle to the ocp.                            |
            %|                                                         |
            %| Return values:                                          |
            %|   n - The number of control inputs.                     |
            %|_________________________________________________________|
            n = numel(obj.controls.mx_vec);
        end
        
        function n = n_p(obj)
            %____________________________________________________________
            %|YOP.OCP/N_P Get the number of free parameters in the OCP. |
            %|                                                          |
            %| Use:                                                     |
            %|   n = n_p(obj)                                           |
            %|   n = obj.n_p()                                          |
            %|   n = obj.n_p                                            |
            %|                                                          |
            %| Parameters:                                              |
            %|   obj - A handle to the ocp.                             |
            %|                                                          |
            %| Return values:                                           |
            %|   n - The number of free parameters.                     |
            %|__________________________________________________________|
            n = numel(obj.parameters.mx_vec);
        end
        
        function n = n_tp(obj)
            %________________________________________________________
            %|YOP.OCP/N_TP Get the number of timepoints in the OCP. |
            %|                                                      |
            %| Use:                                                 |
            %|   n = n_tp(obj)                                      |
            %|   n = obj.n_tp()                                     |
            %|   n = obj.n_tp                                       |
            %|                                                      |
            %| Parameters:                                          |
            %|   obj - A handle to the ocp.                         |
            %|                                                      |
            %| Return values:                                       |
            %|   n - The number of timepoints.                      |
            %|______________________________________________________|
            n = numel(obj.timepoints.mx_vec);
        end
        
        function n = n_int(obj)
            %________________________________________________________
            %|YOP.OCP/N_INT Get the number of integrals in the OCP. |
            %|                                                      |
            %| Use:                                                 |
            %|   n = n_int(obj)                                     |
            %|   n = obj.n_int()                                    |
            %|   n = obj.n_int                                      |
            %|                                                      |
            %| Parameters:                                          |
            %|   obj - A handle to the ocp.                         |
            %|                                                      |
            %| Return values:                                       |
            %|   n - The number of integrals.                       |
            %|______________________________________________________|
            n = numel(obj.integrals.mx_vec);
        end
        
        function [bool, T] = fixed_horizon(obj)
            %______________________________________________________________
            %|YOP.OCP/FIXED_HORIZON Test if the problem horizon is fixed. |
            %|                                                            |
            %| Use:                                                       |
            %|   bool = fixed_horizon(obj)                                |
            %|   bool = obj.fixed_horizon()                               |
            %|   [bool, T] = fixed_horizon(obj)                           |
            %|   [bool, T] = obj.fixed_horizon()                          |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - A handle to the ocp.                               |
            %|                                                            |
            %| Return values:                                             |
            %|   bool - A boolean value (true/false) telling whether      |
            %|          the horizon is fixed.                             |
            %|   T - The fixed horizon length. If the horizon is not      |
            %|       fixed this value should be disregarded.              |
            %|____________________________________________________________|
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
                fprintf(char(propagate_value(obj.differetial_eqs{k}.lhs)));
                fprintf(') == ');
                rhs = propagate_value(obj.differetial_eqs{k}.rhs);
                if length(char(rhs)) > 40
                    vs = yop.get_vars(obj.differetial_eqs{k}.rhs);
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