classdef vnf_data < handle
    properties
        vn = {}
        nv = {}
        vv = {}
        ve = {}
        ev = {}
        ee = {}
        
    end
    methods
        
        function obj = add_vn(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vn = {obj.vn{:}, e{k}};
                end
            else
                obj.vn = {obj.vn{:}, e};
            end
        end
        
        function obj = add_nv(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.nv = {obj.nv{:}, e{k}};
                end
            else
                obj.nv = {obj.nv{:}, e};
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
        
    end
end