classdef srf_lt < yop.srf
    
    properties (Constant)
        name = 'srf_lt'
    end
    
    methods
        function obj = srf_lt(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end