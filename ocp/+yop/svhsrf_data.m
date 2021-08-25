classdef svhsrf_data < handle
    properties
        ve = {}
        ev = {}
    end
    methods
        
        function obj = add_ve(obj, ve)
            obj.ve = {obj.ve{:}, ve};
        end
        
        function obj = add_ev(obj, ev)
            obj.ev = {obj.ev{:}, ev};
        end
        
    end
end