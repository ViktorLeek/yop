classdef ast_eq < yop.ast_relation
    
    properties
        m_ode = false
        m_alg = false
    end
    
    properties (Constant)
        m_name = 'eq'
    end
    
    methods
        function obj = ast_eq(lhs, rhs, ishard, isode, isalg)
            if isa(lhs, 'function_handle')
                val = eq(lhs(1), value(rhs));
                num = eq(lhs(1), numval(rhs));
            elseif isa(rhs, 'function_handle')
                val = eq(value(lhs), rhs(1));
                num = eq(numval(lhs), rhs(1));
            else
                val = eq(value(lhs), value(rhs));
                num = eq(numval(lhs), numval(rhs));
            end
            obj@yop.ast_relation(val, num, lhs, rhs);
            if nargin > 2
                obj.m_hard = ishard;
                obj.m_ode = isode;
                obj.m_alg = isalg;
            end
        end
        
        function obj = ode(obj)
            obj.m_ode = true;
        end
        
        function obj = alg(obj)
            obj.m_alg = true;
        end
        
        function bool = isa_alg(obj)
            bool = obj.m_alg;
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_eq( ...
                lhs, rhs, obj.m_hard, obj.m_ode, obj.m_alg);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_lhs-obj.m_rhs, 0);
        end

    end
    
end