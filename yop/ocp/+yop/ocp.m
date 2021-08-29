classdef ocp < handle
    properties
        % Misc
        name
        
        % Variables
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        
        % odes and algebraic eq's
        alg_eqs
        odes
        
        % objective
        objective
        
        % Unparsed constraints
        constraints
        
        % Parsed constraints
        equality = {}; % Temporary cell
        inequality = {}; % Temporary cell
    end
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
            obj.objective = -objective;  % negate to get min problem
        end
        
        function obj = st(obj, varargin)
            obj.constraints = varargin;
        end
        
        function obj = build(obj)
            % Reset problem before building?
            obj.parse_variables().parse_constraints().set_box_bounds();
        end
        
        function obj = parse_variables(obj)
            vars = get_variables({obj.objective, obj.constraints{:}});
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
        end
        
        function obj = parse_constraints(obj)
            srf = yop.to_srf(obj.constraints);
            hsrf = yop.to_hsrf(srf.get_relations());
            vnf = yop.to_vnf(hsrf);
            dtp = yop.to_dtp(vnf);
            
            % Box constraints: var-num
            for k=1:length(dtp.vn_t)
                bc = dtp.vn_t{k};
                var_expr = bc.lhs;
                bnd = dummy_evaluate(bc.rhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ub(re.idx_var) = bnd(re.idx_expr);
                        var.lb(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: var-num, t==t0
            for k=1:length(dtp.vn_t0)
                bc = dtp.vn_t0{k};
                var_expr = bc.lhs;
                bnd = dummy_evaluate(bc.rhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ub0(re.idx_var) = bnd(re.idx_expr);
                        var.lb0(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub0(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb0(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: var-num, t==tf
            for k=1:length(dtp.vn_tf)
                bc = dtp.vn_tf{k};
                var_expr = bc.lhs;
                bnd = dummy_evaluate(bc.rhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % var == bnd
                        var.ubf(re.idx_var) = bnd(re.idx_expr);
                        var.lbf(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ubf(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lbf(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var
            for k=1:length(dtp.nv_t)
                bc = dtp.nv_t{k};
                var_expr = bc.rhs;
                bnd = dummy_evaluate(bc.lhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ub(re.idx_var) = bnd(re.idx_expr);
                        var.lb(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var, t==t0
            for k=1:length(dtp.nv_t0)
                bc = dtp.nv_t0{k};
                var_expr = bc.rhs;
                bnd = dummy_evaluate(bc.lhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ub0(re.idx_var) = bnd(re.idx_expr);
                        var.lb0(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lb0(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ub0(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Box constraints: num-var, t==tf
            for k=1:length(dtp.nv_tf)
                bc = dtp.nv_tf{k};
                var_expr = bc.rhs;
                bnd = dummy_evaluate(bc.lhs);
                re = yop.reaching_elements(var_expr);
                var = obj.find_variable(re.var.id);
                switch class(bc)
                    case 'yop.ast_eq' % bnd == var
                        var.ubf(re.idx_var) = bnd(re.idx_expr);
                        var.lbf(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_le', 'yop.ast_lt'} % bnd < var
                        var.lbf(re.idx_var) = bnd(re.idx_expr);
                    case {'yop.ast_ge', 'yop.ast_gt'} % bnd > var
                        var.ubf(re.idx_var) = bnd(re.idx_expr);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            
            % Path constraints
            pc = dtp.get_pathcon();
            for k=1:length(pc)
                switch class(pc{k})
                    case 'yop.ast_eq'
                        c = yop.ast_eq(pc{k}.lhs - pc{k}.rhs, 0);
                        obj.add_equality(c);
                    case {'yop.ast_le', 'yop.ast_lt'}
                        c = yop.ast_le(pc{k}.lhs - pc{k}.rhs, 0);
                        obj.add_inequality(c);
                    case {'yop.ast_ge', 'yop.ast_gt'}
                        c = yop.ast_le(pc{k}.rhs - pc{k}.lhs, 0);
                        obj.add_inequality(c);
                end
            end
        end
        
        function obj = set_box_bounds(obj)
            obj.set_box_bnd(...
                'independent', 'lb', 'ub', ...
                yop.defaults().independent_lb, ...
                yop.defaults().independent_ub ...
                );
            
            obj.set_box_bnd(...
                'independent_initial', 'lb', 'ub', ...
                yop.defaults().independent_lb, ...
                yop.defaults().independent_ub ...
                );
            
            obj.set_box_bnd(...
                'independent_final', 'lb', 'ub', ...
                yop.defaults().independent_lb, ...
                yop.defaults().independent_ub ...
                );
            
            obj.set_box_bnd(...
                'states', 'lb', 'ub', ...
                yop.defaults().state_lb, ...
                yop.defaults().state_ub ...
                );
            
            obj.set_box_bnd(...
                'states', 'lb0', 'ub0', ...
                yop.defaults().state_lb0, ...
                yop.defaults().state_ub0 ...
                );
            
            obj.set_box_bnd(...
                'states', 'lbf', 'ubf', ...
                yop.defaults().state_lbf, ...
                yop.defaults().state_ubf ...
                );
            
            obj.set_box_bnd(...
                'algebraics', 'lb', 'ub', ...
                yop.defaults().algebraic_lb, ...
                yop.defaults().algebraic_ub ...
                );
            
            obj.set_box_bnd(...
                'controls', 'lb', 'ub', ...
                yop.defaults().control_lb, ...
                yop.defaults().control_ub ...
                );
            
            obj.set_box_bnd(...
                'controls', 'lb0', 'ub0', ...
                yop.defaults().control_lb0, ...
                yop.defaults().control_ub0 ...
                );
            
            obj.set_box_bnd(...
                'controls', 'lbf', 'ubf', ...
                yop.defaults().control_lbf, ...
                yop.defaults().control_ubf ...
                );
            
            obj.set_box_bnd(...
                'parameters', 'lb', 'ub', ...
                yop.defaults().parameter_lb, ...
                yop.defaults().parameter_ub ...
                );
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
        
        function vars = variables(obj)
            % Find variable might increase in speed depending on the order
            % the variables are conctatenated.
            vars = [...
                obj.independent(:).', ...
                obj.independent_initial(:).', ...
                obj.independent_final(:).', ...
                obj.states(:).', ...
                obj.algebraics(:).', ...
                obj.controls(:).', ...
                obj.parameters(:).' ...
                ];
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
                obj.independent = yop.ocp_variable(t);
            else
                obj.independent(end+1) = yop.ocp_variable(t);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_variable(t);
            else
                obj.independent_initial(end+1) = ...
                    yop.ocp_variable(t);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_variable(t);
            else
                obj.independent_final(end+1) = ...
                    yop.ocp_variable(t);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ocp_variable(x);
            else
                obj.states(end+1) = yop.ocp_variable(x);
            end
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_variable(z);
            else
                obj.algebraics(end+1) = yop.ocp_variable(z);
            end
        end
        
        function obj = add_algebraic_eq(obj, eq)
            if isempty(obj.alg_eqs)
                obj.alg_eqs = yop.ocp_variable(eq);
            else
                obj.alg_eqs(end+1) = yop.ocp_variable(eq);
            end
        end
        
        function obj = add_control(obj, u)
            if isempty(obj.controls)
                obj.controls = yop.ocp_variable(u);
            else
                obj.controls(end+1) = yop.ocp_variable(u);
            end
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_variable(p);
            else
                obj.parameters(end+1) = yop.ocp_variable(p);
            end
        end
        
        function obj = add_inequality(obj, pc)
            obj.inequality = {obj.inequality{:}, pc};
        end
        
        function obj = add_equality(obj, pc)
            obj.equality = {obj.equality{:}, pc};
        end
        
        function present(obj)
            
            for k=obj.variables
                k.store_value();
            end
            
            for t=obj.independent
                t.set_value(sym(t.var.name));
            end
            
            for t0=obj.independent_initial
                t0.set_value(sym(t0.var.name));
            end
            
            for tf=obj.independent_final
                tf.set_value(sym(tf.var.name));
            end
            
            for x=obj.states
                if isscalar(x.var.size)
                    x.set_value(sym(x.var.name));
                else
                    x.set_value(sym(x.var.name, size(x.var)));
                end
            end
            
            for z=obj.algebraics
                if isscalar(z.var.size)
                    z.set_value(sym(z.var.name));
                else
                    z.set_value(sym(z.var.name, size(z.var)));
                end
            end
            
            for u=obj.controls
                if isscalar(u.var.size)
                    u.set_value(sym(u.var.name));
                else
                    u.set_value(sym(u.var.name, size(u.var)));
                end
            end
            
            for p=obj.parameters
                if isscalar(p.var.size)
                    p.set_value(sym(p.var.name));
                else
                    p.set_value(sym(p.var.name, size(p.var)));
                end
            end
            
            % objective function
            x = forward_evaluate(obj.objective);
            if isempty(obj.name)
                title = 'Optimal Control Problem';
            else
                title = obj.name;
            end
            fprintf(['[Yop] ', title, '\n']);
            fprintf('  min\t');
            fprintf(char(x));
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
                fprintf('  Parameters-\n');
            end
            obj.print_box('parameters');
            
            
            if ~isempty(obj.equality)
                fprintf('  Equality\n');
            end
            for k=1:length(obj.equality)
                fprintf('\t');
                fprintf(char(forward_evaluate(obj.equality{k})));
                fprintf('\n');
            end
            fprintf('\n');
            
            if ~isempty(obj.inequality)
                fprintf('  Inequality\n');
            end
            for k=1:length(obj.inequality)
                fprintf('\t');
                fprintf(char(forward_evaluate(obj.inequality{k})));
                fprintf('\n');
            end
            
            for k=obj.variables
                k.restore_value();
            end
            
        end
        
        function obj = print_box(obj, varstr)
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