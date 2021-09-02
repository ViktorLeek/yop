classdef tri_data < handle
    % trancscription invariant data
    properties
        eq_inv = {}  % equality constraint, transcription invariant
        ieq_inv = {} % inequality constraint, transcription invariant
        eq_var = {}  % equality constraint, transcription variant
        ieq_var = {} % inequality constraint, transcription variant
    end
    methods
        
        function obj = add_eq_inv(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.eq_inv = {obj.eq_inv{:}, e{k}};
                end
            else
                obj.eq_inv = {obj.eq_inv{:}, e};
            end
        end
        
        function obj = add_ieq_inv(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ieq_inv = {obj.ieq_inv{:}, e{k}};
                end
            else
                obj.ieq_inv = {obj.ieq_inv{:}, e};
            end
        end        
        
        function obj = add_eq_var(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.eq_var = {obj.eq_var{:}, e{k}};
                end
            else
                obj.eq_var = {obj.eq_var{:}, e};
            end
        end
        
        function obj = add_ieq_var(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ieq_var = {obj.ieq_var{:}, e{k}};
                end
            else
                obj.ieq_var = {obj.ieq_var{:}, e};
            end
        end
        
    end
end