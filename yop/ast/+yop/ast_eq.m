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
            obj@yop.ast_relation(eq(value(lhs), value(rhs)), lhs, rhs);
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