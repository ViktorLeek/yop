classdef binary_expression < yop.node
    
    properties
        lhs
        rhs
    end
    
    methods
        function obj = binary_expression(lhs, rhs)
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
    end
    
    methods % Printing
        
        function xprint(obj, s)
            fprintf([s, '\n']);
            
            begin_child(obj);
            yop.print(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            yop.print(obj.rhs);
            end_child(obj);
        end
        
    end
end