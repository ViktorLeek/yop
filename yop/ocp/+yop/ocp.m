classdef ocp < handle

    properties
        name 
        
        objective
        constraints
        
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        
        ode
        special_nodes
        timepoints
        integrals
        
        equality_constraints
        inequality_constraints
    end
    
    properties (Hidden)
        %         built = false;
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
        
        function [sol, nlp] = solve(obj, varargin)
            ip = inputParser();
            ip.FunctionName = "yop.ocp/solve";
            ip.addParameter('intervals', yop.defaults.control_invervals);
            ip.addParameter('degree', yop.defaults.polynomial_degree);
            ip.addParameter('points', yop.defaults.collocation_points);
            ip.parse(varargin{:});
            
            N  = ip.Results.intervals;
            d  = ip.Results.degree;
            cp = ip.Results.points;
            
            obj.to_canonical();
            
            nlp = yop.direct_collocation(obj, N, d, cp);
            
            solver = casadi.nlpsol( ...
                'solver', ...
                'ipopt', ...
                struct('f', nlp.J, 'x', nlp.w, 'g', [nlp.g; nlp.h]) ...
                );
            
            w0 = ones(size(nlp.w));
            sol = solver( ...
                'x0', w0, ...
                'lbx', nlp.w_lb, ...
                'ubx', nlp.w_ub, ...
                'ubg', [nlp.g_ub; nlp.h_ub], ...
                'lbg', [nlp.g_lb; nlp.h_lb] ...
                );
            
            
            
        end
        
        function obj = to_canonical(obj)
            
            [vars, tps, ints, ders, sn] = yop.ocp.find_special_nodes(...
                obj.objective.ast, obj.constraints{:});
            obj.timepoints = tps;
            obj.integrals = ints;
            obj.special_nodes = sn;
            
            obj.classify_variables(vars);
            
            srf = yop.ocp.to_srf(obj.constraints);
            [box, nbox] = yop.ocp.separate_box(srf);
            obj.set_box(box);
            
            [odes, algs, eqs, ieqs] = yop.ocp.sort_nonbox(nbox);
            obj.vectorize_ode(odes);
            
            obj.equality_constraints = yop.ocp_expr.empty(1,0);
            for k=1:length(eqs)
                obj.equality_constraints(end+1) = ...
                    yop.ocp_expr(eqs{k}.lhs, is_hard(eqs{k}));
            end
            
            obj.inequality_constraints = yop.ocp_expr.empty(1,0);
            for k=1:length(ieqs)
                obj.inequality_constraints(end+1) = ...
                    yop.ocp_expr(ieqs{k}.lhs, is_hard(ieqs{k}));
            end
            
            obj.set_ode_function();
            obj.set_special_functions(tps, ints, ders);
            obj.set_objective_function(tps, ints, ders);
            obj.set_pathcon_functions(tps, ints, ders);
            
            % Error/semantic checking
        end
        
        function obj = vectorize_ode(obj, odes)
            %______________________________________________________________
            %|YOP.OCP/VECTORIZE_ODE Vectorize all ODEs and set default    |
            %|odes for states that are not bound by an ode.               |
            %|                                                            |
            %| Use:                                                       |
            %|   vectorize_ode(obj, odes)                                 |
            %|   obj.vectorize_ode(odes)                                  |
            %|                                                            |
            %| Description:                                               |
            %|   Analyzes which of the state elements the reaches the ode |
            %|   rhs: der(var) == expr. The elements of the entire state  |
            %|   vector that does not reach is given the default          |
            %|   derivative: der(var) == 0.                               |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the ocp.                                 |
            %|   odes - Cell array with canonicalized ODEs.               |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
            ode_var = []; ode_expr = [];
            for k=1:length(odes)
                ode_var = [ode_var(:); odes{k}.lhs(:)];
                ode_expr = [ode_expr(:); odes{k}.rhs(:)];
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
            
            i2e = obj.compute_permutation_vector(ode_var);
            
            ast = yop.ast_eq(ode_var(i2e), ode_expr(i2e));
            obj.ode = yop.ocp_rel(ast);
        end
        
        function [i2e, e2i] = compute_permutation_vector(obj, ode_var)
            %______________________________________________________________!!!!!!UPDATE!!!!!
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
            e2i = output_idx;
            i2e = idx;
        end
        
        function obj = set_ode_function(obj)
            args = { ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters) ...
                };
            set_mx(obj.variables);
            mx_expr = fw_eval(obj.ode.ast.rhs);
            obj.ode.fn = casadi.Function('ode', args, {mx_expr});
        end
        
        function obj = set_special_functions(obj, tps, ints, ders)
            args = { ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx(obj.variables);
            set_mx([tps, ints, ders]);
        
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
        end
        
        function obj = set_objective_function(obj, tps, ints, ders)
            % The objective should be a functino in tps and ints only.
            args = {mx_vec(tps), mx_vec(ints)};
            
            set_mx(obj.variables); % To avoid undefined behaviour.
            set_mx([tps, ints, ders]);
        
            mx_expr = fw_eval(obj.objective.ast);
            obj.objective.fn = casadi.Function('fn', args, {mx_expr});
        end
        
        function obj = set_pathcon_functions(obj, tps, ints, ders)
            
            args = { ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx(obj.variables);
            set_mx([tps, ints, ders]);
            
            for pc = [obj.equality_constraints, obj.inequality_constraints]
                mx_expr = fw_eval(pc.ast);
                pc.fn = casadi.Function('fn', args, {mx_expr});
            end
        end
    end

    %% Variables
    methods
        
        function obj = classify_variables(obj, vars)
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
            %|   vars - Cell array with ast_variable-s.                   |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the ocp.                                 |
            %|____________________________________________________________|
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
            % them, and also to be able to set bounds on the independent
            % variable.
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
        
        function obj = set_box(obj, box)
            [box_t, box_t0, box_tf] = yop.ocp.timed_box(box);
            obj.parse_box(box_t, 'lb', 'ub');
            obj.parse_box(box_t0, 'lb0', 'ub0');
            obj.parse_box(box_tf, 'lbf', 'ubf');
            obj.set_box_bounds(); % Includes processing default values
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
            if isempty(obj.states)
                obj.states = yop.ocp_var.empty(1,0);
            end
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
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_var.empty(1,0);
            end
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
            if isempty(obj.controls)
                obj.controls = yop.ocp_var.empty(1,0);
            end
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
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_var.empty(1,0);
            end
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
    end
    
    %% Internal methods and helper methods
    
    methods
        
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
        
        
        function [vars,tps,ints,ders,sorted] = find_special_nodes(varargin)
            vars = {};
            sorted = yop.ocp_expr.empty(1,0); % Topological order
            tps   = yop.ocp_expr.empty(1,0);
            ints  = yop.ocp_expr.empty(1,0);
            ders  = yop.ocp_expr.empty(1,0);
            visited = [];
            for k=1:length(varargin)
                % Sort all nodes
                [tsort, n_elem, visited] = ...
                    topological_sort(varargin{k}, visited);
                
                % Find special nodes, maintin topological order in sorted.
                for n=1:n_elem
                    if isa(tsort{n}, 'yop.ast_variable')
                        vars{end+1} = tsort{n};
                    elseif isa(tsort{n}, 'yop.ast_timepoint')
                        sn = yop.ocp_expr(tsort{n}, yop.ocp_expr.tp);
                        tps(end+1) = sn;
                        sorted(end+1) = sn;
                    elseif isa(tsort{n}, 'yop.ast_int')
                        sn = yop.ocp_expr(tsort{n}, yop.ocp_expr.int);
                        ints(end+1) = sn;
                        sorted(end+1) = sn;
                    elseif isa(tsort{n}, 'yop.ast_der')
                        %sn = yop.ocp_expr(node, yop.ocp_expr.der);
                        %obj.derivatives(end+1) = sn;
                        %obj.special_nodes(end+1) = sn;
                    end
                end
            end
        end
        
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
            %|YOP.OCP.SORT_NONBOX Sorts the constraints that are not box  |
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
                    eq{end+1} = eq_k; % is canonicalized
                    
                else
                    ieq{end+1} = canonicalize(nbox{k});
                    
                end
            end
            ode = ode(~cellfun('isempty', ode));
            alg = alg(~cellfun('isempty', alg));
            ieq = ieq(~cellfun('isempty', ieq));
            eq  =  eq(~cellfun('isempty', eq));
            
            eq = yop.ocp.split_hard(eq);
            ieq = yop.ocp.split_hard(ieq);
            eq = yop.ocp.split_transcription_invariance(eq);
            ieq = yop.ocp.split_transcription_invariance(ieq);
        end
        
        function con = split_hard(con)
            %______________________________________________________________
            %|YOP.OCP.SPLIT_HARD Split regular and hard constraints.      |
            %|Is only relevant for direct collocation. Hard constraints   |
            %|are applied to every discretization point.                  |
            %|                                                            |
            %| Use:                                                       |
            %|   con = yop.ocp.split_hard(con)                            |
            %|                                                            |
            %| Parameters:                                                |
            %|   con - Cell array with canonlicalized non-box constraints.|
            %|                                                            |
            %| Return values:                                             |
            %|   trcon - Cell array with the constraints separated based  |
            %|           if they are hard or not.                         |
            %|____________________________________________________________|
            hcon = {};
            for k=1:length(con)
                hard = is_hard(con{k});
                if all(hard) || all(~hard)
                    hcon{end+1} = con{k};
                else
                    hcon{end+1} = yop.get_subrel(con{k}, hard);
                    hcon{end+1} = yop.get_subrel(con{k}, ~hard);
                end
            end
            con = hcon;
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
            %|   con - Cell array with canonlicalized non-box constraints.|
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
                    trcon{end+1} = yop.get_subrel(con{k}, invariant);
                    trcon{end+1} = yop.get_subrel(con{k}, ~invariant);
                end
            end
        end
        
    end
    
    
    %% Transcription
    methods
        
        function n = n_x(obj)
            n = n_elem(obj.states);
        end
        
        function n = n_u(obj)
            n = n_elem(obj.controls);
        end
        
        function n = n_p(obj)
            n = n_elem(obj.parameters);
        end
        
        function n = n_tp(obj)
            n = n_elem(obj.timepoints);
        end
        
        function n = n_int(obj)
            n = n_elem(obj.integrals);
        end
                
        function bd = t0_ub(obj)
            bd = obj.independent_initial.ub;
        end
        
        function bd = t0_lb(obj)
            bd = obj.independent_initial.lb;
        end
        
        function bd = tf_ub(obj)
            bd = obj.independent_final.ub;
        end
        
        function bd = tf_lb(obj)
            bd = obj.independent_final.lb;
        end
        
        function bd = x0_ub(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).ub0(:)];
            end
        end
        
        function bd = x0_lb(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).lb0(:)];
            end
        end
        
        function bd = x_ub(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).ub(:)];
            end
        end
        
        function bd = x_lb(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).lb(:)];
            end
        end
        
        function bd = xf_ub(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).ubf(:)];
            end
        end
        
        function bd = xf_lb(obj)
            bd=[];
            for k=1:length(obj.states)
                bd = [bd(:); obj.states(k).lbf(:)];
            end
        end
        
        function bd = u0_ub(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).ub0(:)];
            end
        end
        
        function bd = u0_lb(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).lb0(:)];
            end
        end
        
        function bd = u_ub(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).ub(:)];
            end
        end
        
        function bd = u_lb(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).lb(:)];
            end
        end
        
        function bd = uf_ub(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).ubf(:)];
            end
        end
        
        function bd = uf_lb(obj)
            bd=[];
            for k=1:length(obj.controls)
                bd = [bd(:); obj.controls(k).lbf(:)];
            end
        end
        
        function bd = p_ub(obj)
            bd=[];
            for k=1:length(obj.parameters)
                bd = [bd(:); obj.parameters(k).ub(:)];
            end
        end
        
        function bd = p_lb(obj)
            bd=[];
            for k=1:length(obj.parameters)
                bd = [bd(:); obj.parameters(k).lb(:)];
            end
        end
        
        function bd = z_ub(obj)
            bd=[];
            for k=1:length(obj.algebraics)
                bd = [bd(:); obj.algebraics(k).ub(:)];
            end
        end
        
        function bd = z_lb(obj)
            bd=[];
            for k=1:length(obj.algebraics)
                bd = [bd(:); obj.algebraics(k).lb(:)];
            end
        end
        
        function [bool, t0, tf] = fixed_horizon(obj)
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
            t0 = obj.t0_lb;
            tf = obj.tf_lb;
        end
        
    end
    
    %% Presentation
    methods
        
        function obj = present(obj)
            obj.build();
            
            obj.independent.set_sym();
            obj.independent_initial.set_sym();
            obj.independent_final.set_sym();
            obj.states.set_sym();
            obj.algebraics.set_sym();
            obj.controls.set_sym();
            obj.parameters.set_sym();
            
            % objective function
            val = propagate_value(obj.objective.ast);
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
            
            
            if ~isempty(obj.ode)
                fprintf('  ODE\n');
            end
            
            fprintf('\tder(');
            fprintf(char(propagate_value(obj.ode.lhs)));
            fprintf(') == ');
            rhs = propagate_value(obj.ode.rhs);
            if length(char(rhs)) > 40
                vs = yop.get_vars(obj.ode.rhs);
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
            
            if ~isempty(obj.equality_constraints)
                fprintf('  Equality\n');
            end
            for k=1:length(obj.equality_constraints)
                fprintf('\t');
                %                 fprintf(char(propagate_value(obj.eq{k})));
                fprintf('\n');
            end
            
            if ~isempty(obj.inequality_constraints)
                fprintf('  Inequality\n');
            end
            for k=1:length(obj.inequality_constraints)
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

% function obj = discretize(obj, tro)
% [sn, n_tp, n_int, n_der] = mx_functionx( ...
%     obj.objective, ...
%     obj.independent, ...
%     obj.states, ...
%     obj.controls, ...
%     obj.parameters ...
%     );
% tro.add_objective(obj.objective, sn, n_tp, n_int, n_der);
%
% mx_function( ...
%     obj.ode.rhs, ...
%     obj.independent, ...
%     obj.states, ...
%     obj.controls, ...
%     obj.parameters ...
%     );
% tro.add_ode(obj.ode.rhs.fn);
%
% for eq = obj.equality_constraints
%     [sn, n_tp, n_int, n_der] = ...
%         mx_functionx( ...
%         eq, ...
%         obj.independent, ...
%         obj.states, ...
%         obj.controls, ...
%         obj.parameters ...
%         );
%     tro.add_eq(eq, sn, n_tp, n_int, n_der);
% end
%
% for ieq = obj.inequality_constraints
%     [sn, n_tp, n_int, n_der] = ...
%         mx_functionx( ...
%         ieq, ...
%         obj.independent, ...
%         obj.states, ...
%         obj.controls, ...
%         obj.parameters ...
%         );
%     tro.add_ieq(ieq, sn, n_tp, n_int, n_der);
% end
% end