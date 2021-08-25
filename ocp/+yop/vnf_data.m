classdef vnf_data < handle
    properties
        vn = {}
        nv = {}
        ve = {}
        ev = {}
    end
    methods
        
        function obj = add_vn(obj, vn)
            obj.vn = {obj.vn{:}, vn};
        end
        
        function obj = add_nv(obj, nv)
            obj.nv = {obj.nv{:}, nv};
        end
        
        function obj = add_ve(obj, ve)
            obj.ve = {obj.ve{:}, ve};
        end
        
        function obj = add_ev(obj, ev)
            obj.ev = {obj.ev{:}, ev};
        end
        
    end
end