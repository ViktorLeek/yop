classdef srf_ne < yop.srf
    
    properties (Constant)
        name = 'srf_ne'
    end
    
    methods
        function obj = srf_ne(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end