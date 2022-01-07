classdef ast_independent < yop.ast_variable
    
    properties (Constant)
        m_reducible = false;
    end
    
    methods
        function obj = ast_independent(name, w, os)
            obj@yop.ast_variable(name, w, os, false, yop.var_type.time);
        end
    end
end