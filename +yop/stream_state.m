classdef stream_state < handle
    
    properties
        indent_level = 1;
        branches = false(1e1,1);
    end
    
    methods
        
        function obj = stream_state()
        end
        
        function reset(obj)
            obj.indent_level = 1;
            obj.branches = false(1e1,1);
        end
        
    end
end