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
        
        function obj = add_vu(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vu = {obj.vu{:}, e{k}};
                end
            else
                obj.vu = {obj.vu{:}, e};
            end
        end
        
        function obj = add_eu(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.eu = {obj.eu{:}, e{k}};
                end
            else
                obj.eu = {obj.eu{:}, e};
            end
        end
        
        function obj = add_vv(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vv = {obj.vv{:}, e{k}};
                end
            else
                obj.vv = {obj.vv{:}, e};
            end
        end
        
        function obj = add_ve(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ve = {obj.ve{:}, e{k}};
                end
            else
                obj.ve = {obj.ve{:}, e};
            end
        end
        
        function obj = add_ev(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ev = {obj.ev{:}, e{k}};
                end
            else
                obj.ev = {obj.ev{:}, e};
            end
        end
        
        function obj = add_ee(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ee = {obj.ee{:}, e{k}};
                end
            else
                obj.ee = {obj.ee{:}, e};
            end
        end
        
        function obj = clear_unknown(obj)
            obj.vu = {};
            obj.eu = {};
        end
    end
end