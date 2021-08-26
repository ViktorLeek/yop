classdef ocp < handle
    properties
        objective
        constraints
        variables
        box
        equality
        inequality
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        alg_eqs  % algebraic equations
        odes
    end
    methods
        function obj = ocp()
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
            obj.parse_variables();
            obj.parse_constraints();
        end
        
        function obj = parse_variables(obj)
            vars = yop.get_variables({obj.objective, obj.constraints{:}});
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
            
            % Box constraints
            for k=1:length(dtp.vn_t)
                obj.add_box(dtp.vn_t{k});
            end
            
            for k=1:length(dtp.vn_t0)
                obj.add_box0(dtp.vn_t0{k});
            end
            
            for k=1:length(dtp.vn_tf)
                obj.add_boxf(dtp.vn_tf{k});
            end
            
            for k=1:length(dtp.nv_t)
                obj.add_box(dtp.nv_t{k});
            end
            
            for k=1:length(dtp.nv_t0)
                obj.add_box0(dtp.nv_t0{k});
            end
            
            for k=1:length(dtp.nv_tf)
                obj.add_boxf(dtp.nv_tf{k});
            end
            
            % Path constraints
        end
        
        function obj = match_box_constraints(obj, vnf)
            for k=1:length(vnf.vn)
                bk = vnf.vn{k};
                obj.match_box_constraint(...
                    class(bk), ...
                    bk.lhs, ...
                    dummy_evaluate(bk.rhs) ...
                    );
            end
            
            % The same for nv
        end
        
%         function obj = match_box_constraint(obj, class, variable, bound)
%             switch class(bk)
%                 case 'yop.ast_eq'
%                     obj.set_box_equality(variable, value, bound);
%                     
%                 case {'yop.ast_gt', 'yop.ast_ge'}
%                     
%                 case {'yop.ast_lt', 'yop.ast_le'}
%                     
%                 otherwise
%                     error(['[Yop] Error: Illegal relation for a ', ...
%                         'box constraint.']);
%             end
%         end
%         
%         function obj = set_box_equality(obj, var_expr, value)
%             % 1) Find the variable which the constraint concerns
%             % 2) Get the reaching indices
%             % 3) Get timepoints
%             % 4) For every element
%             %    1) Om det inte är en tp
%             %         - Sätt övre och undre gräns till value
%             %    2) Om det är en tp
%             %         - Sätt övre och undre gräns för tidpunkten
%             %
%             
%             re = reaching_elements(var_expr);
%             assert(length(re)==1, '[Yop] Unexpected error.');
%             
%             var = obj.find_variable(re.var.id);
%             [is_tp, tps] = isa_timepoint(var_expr);
%             
%             for k=1:length(is_tp)
%                 if is_tp(k)
%                     switch tps(k)
%                         case yop.initial_timepoint()
%                             % try-catch? ub0 may not exist
%                             var.ub0(re.idx_var) = value(re.idx_expr);
%                             var.lb0(re.idx_var) = value(re.idx_expr);
%                             
%                         case yop.final_timepoint()
%                             var.ubf(re.idx_var) = value(re.idx_expr);
%                             var.lbf(re.idx_var) = value(re.idx_expr);
%                             
%                         otherwise
%                             % Här har det smugit sig in ett inequality
%                             % constraint.
%                             disp('ieq')
%                     end
%                 else
%                     var.ub(re.idx_var) = value(re.idx_expr);
%                     var.lb(re.idx_var) = value(re.idx_expr);
%                     
%                 end
%             end
%             
%         end
        
        function var = find_variable(obj, id)
            for k=1:length(obj.independent)
                if obj.independent{k}.id == id
                    var = obj.independent{k};
                    return
                end
            end
            
            for k=1:length(obj.independent_initial)
                if obj.independent_initial{k}.id == id
                    var = obj.independent_initial{k};
                    return
                end
            end
            
            for k=1:length(obj.independent_final)
                if obj.independent_final{k}.id == id
                    var = obj.independent_final{k};
                    return
                end
            end
            
            for k=1:length(obj.states)
                if obj.states{k}.id == id
                    var = obj.states{k};
                    return
                end
            end
            
            for k=1:length(obj.algebraics)
                if obj.algebraics{k}.id == id
                    var = obj.algebraics{k};
                    return
                end
            end
            
            for k=1:length(obj.controls)
                if obj.controls{k}.id == id
                    var = obj.controls{k};
                    return
                end
            end
            
            for k=1:length(obj.parameters)
                if obj.parameters{k}.id == id
                    var = obj.parameters{k};
                    return
                end
            end
            
            error('[Yop] Error: ID not found');
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ocp_independent(t);
            else
                obj.independent(end+1) = yop.ocp_independent(t);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_independent_initial(t);
            else
                obj.independent_initial(end+1) = ...
                    yop.ocp_independent_initial(t);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_independent_final(t);
            else
                obj.independent_final(end+1) = ...
                    yop.ocp_independent_final(t);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ocp_state(x);
            else
                obj.states(end+1) = yop.ocp_state(x);
            end
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_algebraic(z);
            else
                obj.algebraics(end+1) = yop.ocp_algebraic(z);
            end
        end
        
        function obj = add_algebraic_eq(obj, eq)
            if isempty(obj.alg_eqs)
                obj.alg_eqs = yop.ocp_algebraic_eq(eq);
            else
                obj.alg_eqs(end+1) = yop.ocp_algebraic_eq(eq);
            end
        end
        
        function obj = add_control(obj, u)
            if isempty(obj.controls)
                obj.controls = yop.ocp_control(u);
            else
                obj.controls(end+1) = yop.ocp_control(u);
            end
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_parameter(p);
            else
                obj.parameters(end+1) = yop.ocp_parameter(p);
            end
        end
        
%         function obj = add_box(obj, bc)
%             obj.box = [obj.box(:)', {bc}];
%         end
%         
%         function obj = add_equality(obj, bc)
%             obj.equality = [obj.equality(:)', {bc}];
%         end
%         
%         function obj = add_inequality(obj, bc)
%             obj.inequality = [obj.inequality(:)', {bc}];
%         end
%         
%         function obj = add_ode(obj, ode)
%             if isempty(obj.odes)
%                 obj.odes = yop.ocp_ode(ode.var, ode.expr);
%             else
%                 obj.odes(end+1) = yop.ocp_ode(ode.var, ode.expr);
%             end
%         end
        
    end
end