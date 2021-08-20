classdef ocp_algebraic_eq < yop.node
    properties
        alg
    end
    methods
        function obj = ocp_algebraic(alg)
            obj@yop.node();
            obj.alg = alg;
        end
        
        function draw(obj)
            fprintf('ocp_algebraic_eq(alg)\n');
            
            last_child(obj);
            draw(obj.alg);
            end_child(obj);
        end
    end
end