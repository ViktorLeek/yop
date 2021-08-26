classdef dtp_data < handle
    properties
        vn_t  = {}
        vn_t0 = {}
        vn_tf = {}
        nv_t  = {}
        nv_t0 = {}
        nv_tf = {}
        vv    = {}
        ve    = {}
        ev    = {}
        ee    = {}
    end
    methods
        
        function obj = add_vn_t(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vn_t = {obj.vn_t{:}, e{k}};
                end
            else
                obj.vn_t = {obj.vn_t{:}, e};
            end
        end
        
        function obj = add_vn_t0(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vn_t0 = {obj.vn_t0{:}, e{k}};
                end
            else
                obj.vn_t0 = {obj.vn_t0{:}, e};
            end
        end
        
        function obj = add_vn_tf(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.vn_tf = {obj.vn_tf{:}, e{k}};
                end
            else
                obj.vn_tf = {obj.vn_tf{:}, e};
            end
        end
        
        function obj = add_nv_t(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.nv_t = {obj.nv_t{:}, e{k}};
                end
            else
                obj.nv_t = {obj.nv_t{:}, e};
            end
        end
        
        function obj = add_nv_t0(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.nv_t0 = {obj.nv_t0{:}, e{k}};
                end
            else
                obj.nv_t0 = {obj.nv_t0{:}, e};
            end
        end
        
        function obj = add_nv_tf(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.nv_tf = {obj.nv_tf{:}, e{k}};
                end
            else
                obj.nv_tf = {obj.nv_tf{:}, e};
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