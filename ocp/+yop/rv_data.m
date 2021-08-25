classdef rv_data < handle
    % Reaching variables data
    properties
        var  % The variable the result concerns
        enum % The enumeration of the elements of the variable
        reaching  % The enumerations that reaches the end expression
        index % The indices in the expression that corresponds to the enum
        tmp_value
    end
    methods
        function idx0 = enumerate(obj, idx0)
            n = prod(size(obj.var));
            obj.enum = idx0:(idx0+n-1);
            idx0 = idx0 + n;
            obj.set_value();
        end
        
        function obj = set_value(obj)
            % Set the value of the underlying variable to the enumeration
            % with the same shape as the variable in order for the
            % evaluation to work as it should.
            obj.tmp_value = obj.var.m_value;
            obj.var.m_value = reshape(obj.enum, size(obj.var));
        end
        
        function obj = compute_indices(obj, re)
            obj.set_reaching(re);
            obj.index = false(size(re));
            for k=1:length(obj.reaching)
                obj.index(re == obj.reaching(k)) = true;
            end
        end
        
        function obj = set_reaching(obj, re)
            % re - elements that reaches the expression
            obj.reaching = re(re >= obj.enum(1) & re <= obj.enum(end));
        end
        
        function obj = restore_value(obj)
            obj.var.m_value = obj.tmp_value;
        end
        
    end
end