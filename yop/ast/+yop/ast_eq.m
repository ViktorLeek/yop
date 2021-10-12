classdef ast_eq < yop.ast_relation
    
    properties
         m_alg = false
    end
    
    properties (Constant)
        name = 'eq'
    end
    
    methods
        function obj = ast_eq(lhs, rhs, ishard, isalg)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(eq(ones(size(lhs)), ones(size(rhs))));
            if nargin > 2
                obj.m_hard = ishard;
                obj.m_alg = isalg;
            end
        end
        
        function value = evaluate(obj)
            value = eq(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = eq(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
        
        function obj = alg(obj)
            obj.m_alg = true;
        end
        
        function bool = is_alg(obj)
            bool = obj.m_alg;
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_eq(lhs, rhs, obj.m_hard, obj.m_alg);
        end
        
        function [isode, ode, non_ode, der_id] = isa_ode(obj)
            der_id = [];
            
            %ISA_ODE
            [isa_der_lhs, ID_lhs] = isa_der(obj.lhs);
            [isa_der_rhs, ID_rhs] = isa_der(obj.rhs);
            isode_lhs = isa_state(obj.lhs) & isa_der_lhs;
            isode_rhs = isa_state(obj.rhs) & isa_der_rhs;
            isode = bitxor(isode_lhs, isode_rhs);
            
            if all(~isode)
                ode = [];
            else
                ode = canonicalize_ode( ...
                    yop.get_subrel(obj, isode), ...
                    isode_lhs(isode), ...
                    isode_rhs(isode) ...
                    );
                idl = ID_lhs(isode);
                idr = ID_rhs(isode);
                der_id = [idl(:); idr(:)]; 
                der_id = der_id(der_id~=0);
            end
            
            if all(isode)
                non_ode = [];
            else
                non_ode = canonicalize(yop.get_subrel(obj, ~isode));
            end
            
            assert( ... that both sides are not ode's
                all((isode_lhs & isode_rhs &~isode)==false),...
                ['[Yop] Error: Cannot have ode on both ', ...
                'sides of ==.']);
        end
        
        function code = canonicalize_ode(ode, is_ode_lhs, is_ode_rhs)
            if all(is_ode_lhs)
                code = ode;
            elseif all(is_ode_rhs)
                fn = get_constructor(ode);
                code = fn(ode.rhs, ode.lhs);
            else
                error('[Yop] Unexpected error.');
            end
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.lhs-obj.rhs, 0);
        end
        
        function cbox = canonicalize_box(box)
            isvar = isa_variable(box.lhs);
            if all(isvar) % var == num
                cbox = box;
            else % num == var
                fn = get_constructor(box);
                cbox = fn(box.rhs, box.lhs);
            end
        end
    end
    
end