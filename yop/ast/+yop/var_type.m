classdef var_type
    properties (Constant) % OCP
        not_var = 0
        
        variables_start = 1
        time = 1
        time0 = 2
        timef = 3
        state = 4
        algebraic = 5
        control = 6
        parameter = 7
        variables_stop = 7
    end
    methods (Static)
        function bool = isa_variable(type)
            bool = ...
                type >= yop.var_type.variables_start && ...
                type <= yop.var_type.variables_stop;
        end
    end
end