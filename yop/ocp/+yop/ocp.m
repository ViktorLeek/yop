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
            
            % 5) convert to distinct timepoint form (final form for boxcon)
            dtp = yop.to_dtp(vnf);
            obj.parse_box_constraints(dtp);
            obj.set_box_bounds();
            
            % 6) Sort the nonbox constraints: ode, inequality, equality
            nbc = yop.sort_nonbox(dtp.get_pathcon());
            obj.odes = nbc.odes;
            
            % 7) Convert to transcription invariant/variant form
            tri = yop.tr_invar(nbc); 
            
            % !!!!!!! Spara in tri information !!!!!!!!!!
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
                bnd = yop.prop_num(bc.rhs);
                re = yop.reaching_elems(var_expr);
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
                bnd = yop.prop_num(bc.rhs);
                re = yop.reaching_elems(var_expr);
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
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
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
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
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
                bnd = yop.prop_num(bc.lhs);
                re = yop.reaching_elems(var_expr);
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
        
%         function obj = add_inequality(obj, pc)
%             obj.inequality = {obj.inequality{:}, pc};
%         end
%         
%         function obj = add_equality(obj, pc)
%             obj.equality = {obj.equality{:}, pc};
%         end
        
    end
    
    %% Transcription
    methods
        
        function n = nx(obj)
            n = 0;
            for k=1:length(obj.states)
                n = n + prod(size(obj.states(k).var));
            end
        end
        
        function n = nu(obj)
            n = 0;
            for k=1:length(obj.controls)
                n = n + prod(size(obj.controls(k).var));
            end
        end
        
        function np = np(obj)
            np = 0;
            for k=1:length(obj.parameters)
                np = np + prod(size(obj.parameters(k).var));
            end
        end
        
        
        function [bool, T] = fixed_horizon(obj)
            bool = true;
            for k=1:length(obj.independent_initial)
                bool = bool && (obj.independent_initial(k).lb == ...
                    obj.independent_initial(k).ub);
            end
            for k=1:length(obj.independent_final)
                bool = bool && (obj.independent_final(k).lb == ...
                    obj.independent_final(k).ub);
            end
            T = obj.independent_final(1).lb-obj.independent_initial(1).lb;
        end
        
        function t_vec = t0(obj)
            t_vec = [];
            for v=obj.independent_initial
                t_vec = [t_vec(:); v.sym(:).'];
            end
        end
        
        function t_vec = tf(obj)
            t_vec = [];
            for v=obj.independent_final
                t_vec = [t_vec(:); v.sym(:).'];
            end
        end
        
        function t_vec = t(obj)
            t_vec = [];
            for v=obj.independent
                t_vec = [t_vec(:); v.sym(:).'];
            end
        end
        
        function x_vec = x(obj)
            x_vec = [];
            for v=obj.states
                x_vec = [x_vec(:); v.sym(:).'];
            end
        end
        
        function u_vec = u(obj)
            u_vec = [];
            for v=obj.controls
                u_vec = [u_vec(:); v.sym(:).'];
            end
        end
        
        function p_vec = p(obj)
            p_vec = [];
            for v=obj.parameters
                p_vec = [p_vec(:); v.sym(:).'];
            end
        end
        
        function fn = expr_fn(obj, expr)
            obj.set_sym();
            e = fw_eval(expr);
            obj.reset_variables();
            fn = casadi.Function('fn', {obj.t, obj.x, obj.u, obj.p}, {e});
        end
        
        function [ode_expr, ode_var] = ode(obj)
            obj.set_sym();
            var = []; expr = [];
            for ode = ocp.odes
                % .var expected short, .expr expected longer
                lhs = evaluate(ode.var);
                rhs = fw_eval(ode.expr);
                var = [var(:); lhs(:).'];
                expr = [expr(:); rhs(:).'];
            end
            obj.reset_sym();
            
            ode_var = casadi.Function('ode_var', ...
                {obj.t, obj.x, obj.u, obj.p}, {var});
            
            ode_expr = casadi.Function('ode_expr', ...
                {obj.t, obj.x, obj.u, obj.p}, {expr});
        end
        
        function [t0_lb, t0_ub, tf_lb, tf_ub] = t_bd(obj)
            t0_lb=[]; t0_ub=[]; tf_lb=[]; tf_ub=[];
            for k=1:length(obj.independent_initial)
                t0_lb = [t0_lb(:); obj.independent_initial(k).lb(:).'];
                t0_ub = [t0_ub(:); obj.independent_initial(k).ub(:).'];
            end
            for k=1:length(obj.independent_final)
                tf_lb = [tf_lb(:); obj.independent_final(k).lb(:).'];
                tf_ub = [tf_ub(:); obj.independent_final(k).ub(:).'];
            end
        end
        
        function [x0_lb, x0_ub, x_lb, x_ub, xf_lb, xf_ub] = x_bd(obj)
            x0_lb=[]; x0_ub=[]; x_lb=[]; x_ub=[]; xf_lb=[]; xf_ub=[];
            for k=1:length(obj.states)
                x_lb = [x_lb(:); obj.states(k).lb(:).'];
                x_ub = [x_ub(:); obj.states(k).ub(:).'];
                x0_lb = [x0_lb(:); obj.states(k).lb0(:).'];
                x0_ub = [x0_ub(:); obj.states(k).ub0(:).'];
                xf_lb = [xf_lb(:); obj.states(k).lbf(:).'];
                xf_ub = [xf_ub(:); obj.states(k).ubf(:).'];
            end
        end
        
        function [u0_lb, u0_ub, u_lb, u_ub, uf_lb, uf_ub] = u_bd(obj)
            u0_lb=[]; u0_ub=[]; u_lb=[]; u_ub=[]; uf_lb=[]; uf_ub=[];
            for k=1:length(obj.controls)
                u_lb = [u_lb(:); obj.controls(k).lb(:).'];
                u_ub = [u_ub(:); obj.controls(k).ub(:).'];
                u0_lb = [u0_lb(:); obj.controls(k).lb0(:).'];
                u0_ub = [u0_ub(:); obj.controls(k).ub0(:).'];
                uf_lb = [uf_lb(:); obj.controls(k).lbf(:).'];
                uf_ub = [uf_ub(:); obj.controls(k).ubf(:).'];
            end
        end
        
        function [p_lb, p_ub] = p_bd(obj)
            p_lb=[]; p_ub=[];
            for k=1:length(obj.parameters)
                p_lb = [p_lb(:); obj.parameters(k).lb(:).'];
                p_ub = [p_ub(:); obj.parameters(k).ub(:).'];
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
            x = fw_eval(obj.objective);
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
                fprintf('  Parameters\n');
            end
            obj.print_box('parameters');
            
            
            if ~isempty(obj.odes)
                fprintf('  ODE\n');
            end
            for k=1:length(obj.odes)
                fprintf('\tder(');
                fprintf(char(fw_eval(obj.odes(k).var)));
                fprintf(') == ');
                rhs = char(fw_eval(obj.odes(k).expr));
                if length(rhs) > 40
                    vs = get_vars(obj.odes(k).expr);
                    args = '';
                    for n=1:length(vs)
                        args = [args(:).', vs{n}.name, ', '];
                    end
                    fprintf(['f(', args(1:end-2), ')']);
                else
                    fprintf(rhs);
                end
                fprintf('\n');
            end
            
            if ~isempty(obj.equality)
                fprintf('  Equality\n');
            end
            for k=1:length(obj.equality)
                fprintf('\t');
                fprintf(char(fw_eval(obj.equality{k})));
                fprintf('\n');
            end
            
            if ~isempty(obj.inequality)
                fprintf('  Inequality\n');
            end
            for k=1:length(obj.inequality)
                fprintf('\t');
                fprintf(char(fw_eval(obj.inequality{k})));
                fprintf('\n');
            end
            
            for k=obj.variables
                k.reset_value();
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


% function obj = set_variables(obj, t, t0, tf, x, u, p)
%             obj.set_independent(t);
%             obj.set_independent_initial(t0);
%             obj.set_independent_final(tf);
%             obj.set_states(x);
%             obj.set_controls(u);
%             obj.set_parameters(p);
%         end
%         
%         function obj = reset_vars(obj)
%             obj.set_independent();
%             obj.set_independent_initial();
%             obj.set_independent_final();
%             obj.set_states();
%             obj.set_controls();
%             obj.set_parameters();
%         end
%         
%         function obj = set_independent(obj, t)
%             for k=1:length(obj.independent)
%                 obj.independent(k).set_value(t);
%             end
%         end
%         
%         function obj = reset_independent(obj)
%             for k=1:length(obj.independent)
%                 obj.independent(k).reset_value();
%             end
%         end
%         
%         function obj = set_independent_initial(obj, t0)
%             for k=1:length(obj.independent_initial)
%                 obj.independent_initial(k).set_value(t0);
%             end
%         end
%         
%         function obj = reset_independent_initial(obj)
%             for k=1:length(obj.independent_initial)
%                 obj.independent_initial(k).reset_value();
%             end
%         end
%         
%         function obj = set_independent_final(obj, tf)
%             for k=1:length(obj.independent_final)
%                 obj.independent_final(k).set_value(tf);
%             end
%         end
%         
%         function obj = reset_independent_final(obj)
%             for k=1:length(obj.independent_final)
%                 obj.independent_final(k).reset_value();
%             end
%         end
%         
%         function obj = set_states(obj, x)
%             for k=1:length(obj.states)
%                 obj.states(k).set_value(x{k});
%             end
%         end
%         
%         function obj = reset_states(obj)
%             for k=1:length(obj.states)
%                 obj.states(k).reset_value();
%             end
%         end
%         
%         function obj = set_algebraics(obj, z)
%             for k=1:length(obj.algebraics)
%                 obj.algebraics(k).set_value(z{k});
%             end
%         end
%         
%         function obj = reset_algebraics(obj)
%             for k=1:length(obj.algebraics)
%                 obj.algebraics(k).set_value();
%             end
%         end
%         
%         function obj = set_controls(obj, u)
%             for k=1:length(obj.controls)
%                 obj.controls(k).set_value(u{k});
%             end
%         end
%         
%         function obj = reset_controls(obj)
%             for k=1:length(obj.controls)
%                 obj.controls(k).reset_value();
%             end
%         end
%         
%         function obj = set_parameters(obj, p)
%             for k=1:length(obj.parameters)
%                 obj.parameters(k).set_value(p{k});
%             end
%         end
%         
%         function obj = reset_parameters(obj)
%             for k=1:length(obj.parameters)
%                 obj.parameters(k).reset_value();
%             end
%         end


% % Path constraints
% pc = dtp.get_pathcon();
% for k=1:length(pc)
%     pck = pc{k};
%     switch class(pck)
%         case 'yop.ast_eq'
%             
%             dl = isa_der(pck.lhs);
%             if all(dl) % der(v) == expr
%                 obj.add_ode(pck.lhs, pck.rhs);
%                 continue;
%             elseif all(~dl) % expr == ??
%                 r = pck; % remaining relations
%             else % (der(x), expr) == (expr, ??)
%                 tmp = yop.get_subrelation(pck, dl);
%                 obj.add_ode(tmp.lhs, tmp.rhs);
%                 r = yop.ast_eq(yop.get_subrelation(pck, ~dl));
%                 % obj.add_ode(pck.lhs(dl), pck.rhs(dl));
%                 %r = yop.ast_eq(pck.lhs(~dl), pck.rhs(~dl));
%             end
%             
%             dr = isa_der(r.rhs);
%             if all(dr)
%                 obj.add_ode(r.rhs, r.lhs);
%             elseif all(~dr)
%                 obj.add_equality(yop.ast_eq(r.lhs-r.rhs, 0));
%             else
%                 tmp = yop.get_subrelation(r, dr);
%                 obj.add_ode(tmp.rhs, tmp.lhs);
%                 obj.add_equality(yop.get_subrelation(r, ~dr));
%                 %obj.add_ode(r.rhs(dr), r.lhs(dr));
%                 %obj.add_equality(...
%                 %    yop.ast_eq(r.lhs(~dr)-r.rhs(~dr), 0));
%             end
%             
%         case {'yop.ast_le', 'yop.ast_lt'}
%             c = yop.ast_le(pck.lhs - pck.rhs, 0);
%             obj.add_inequality(c);
%         case {'yop.ast_ge', 'yop.ast_gt'}
%             c = yop.ast_le(pck.rhs - pck.lhs, 0);
%             obj.add_inequality(c);
%     end
% end