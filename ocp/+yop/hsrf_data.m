classdef hsrf_data < handle
    properties
        vu = {} % variables left, unknown right
        eu = {} % expressions left, unknown right
        vv = {}
        ve = {}
        ev = {}
        ee = {}
    end
    methods
        function obj = add_vu(obj, vu)
            obj.vu = {obj.vu{:}, vu};
        end
        
        function obj = add_eu(obj, eu)
            obj.eu = {obj.eu{:}, eu};
        end
        
        function obj = add_vv(obj, vv)
            obj.vv = {obj.vv{:}, vv};
        end
        
        function obj = add_ve(obj, ve)
            obj.ve = {obj.ve{:}, ve};
        end
        
        function obj = add_ev(obj, ev)
            obj.ev = {obj.ev{:}, ev};
        end
        
        function obj = add_ee(obj, ee)
            obj.ee = {obj.ee{:}, ee};
        end
        
        function obj = clear_unknown(obj)
            obj.vu = {};
            obj.eu = {};
        end
    end
end