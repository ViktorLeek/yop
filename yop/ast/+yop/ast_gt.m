classdef ast_gt < yop.ast_relation
    
    properties (Constant)
        name = 'gt'
    end
    
    methods
        function obj = ast_gt(lhs, rhs, ishard)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(gt(ones(size(lhs)), ones(size(rhs))));
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function value = evaluate(obj)
            value = gt(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = gt(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_gt(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.rhs-obj.lhs, 0);
        end
        
        function cbox = canonicalize_box(box)
            isvar = isa_variable(box.lhs);
            if all(isvar) % var >= num
                cbox = box;
            else % num >= var
                cbox = yop.ast_le(box.rhs, box.lhs, box.m_hard);
            end
        end
    end
    
end