classdef srf_eq < yop.srf
    
    properties (Constant)
        name = 'srf_eq'
    end
    
    methods
        function obj = srf_eq(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
        
        function [type, constraint] = get_constraint(obj)
            if isa_variable(obj.lhs) && isnumeric(obj.rhs)
                type = 'box';
                constraint = yop.box_equality(obj.lhs, obj.rhs);
                
            elseif isnumeric(obj.lhs) && isa_variable(obj.rhs)
                type = 'box';
                constraint = yop.box_equality(obj.rhs, obj.lhs);
                
            elseif is_differential(obj.lhs)
                type = 'differential';
                constraint = yop.differential(obj.lhs, obj.rhs);
                
            elseif is_differential(obj.rhs)
                type = 'differential';
                constraint = yop.differential(obj.rhs, obj.lhs);
                
            elseif is_algebraic(obj.lhs)
                type = 'algebraic';
                constraint = yop.algebraic(obj.lhs-obj.rhs);
                
            elseif is_algebraic(obj.rhs)
                type = 'algebraic';
                constraint = yop.algebraic(obj.lhs-obj.rhs);
                
            else
                type = 'equality';
                constraint = yop.equality(obj.lhs-obj.rhs);
            end
        end
    end
end