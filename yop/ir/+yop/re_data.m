classdef re_data < handle
    % Reaching variables data
    properties
        var  % The variable the result concerns
        enum % The enumeration of the elements of the variable
        reaching  % The enumerations that reaches the end expression
        idx_expr % The indices in the expression that corresponds to the enum
    end
    methods
        function [idx0, idx] = enumerate(obj, idx0)
            % Enumerate all elements, starting with idx0
            n = prod(size(obj.var));
            idx = idx0:(idx0+n-1);
            obj.enum = idx;
            obj.set_value();
            idx0 = idx0 + n;
        end
        
        function obj = set_value(obj)
            % Set the value of the underlying variable to the enumeration
            % with the same shape as the variable in order for the
            % evaluation to work as it should.
            obj.var.m_value = reshape(obj.enum, size(obj.var));
        end
        
        function obj = set_expr_idx(obj, re)
            % Computes the index in the expression the variable reaches.
            % It returns a logical array with one element for every element
            % of the exprssion. 'expr(re_data(k).idx_expr)' gives the
            % elements of variable k that reaches the expression.
            obj.set_reaching(re);
            obj.idx_expr = false(size(re));
            for k=1:length(obj.reaching)
                obj.idx_expr(re == obj.reaching(k)) = true;
            end
        end
        
        function idx = idx_var(obj)
            % The indices of the variable that reaches the expression. 
            % This means that re_data(k).var(re_data(k).idx_var) reaches 
            % the expression.
            idx = obj.reaching - obj.enum(1) + 1;
        end
        
        function obj = set_reaching(obj, re)
            % re - elements that reaches the expression
            obj.reaching = re(re >= obj.enum(1) & re <= obj.enum(end));
        end
        
    end
end