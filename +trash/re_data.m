classdef re_data < handle
    %______________________________________________________________________
    %|YOP.RE_DATA Helper data type for doing reaching elements analysis.  |
    %|                                                                    |
    %| Use:                                                               |
    %|   re_data = yop.re_data();                                         |
    %|                                                                    |
    %| Description:                                                       |
    %|   The idea is to enumerate all elements and then evaluate the      |
    %|   expression with the enumeration as input. Since it is possible   |
    %|   to test which elements of the expression that corresponds to     |
    %|   a variable without this analysis, it is possible to remove the   |
    %|   elements that are not enumerations, and by that know how the     |
    %|   the elements are permuted by the expression. The type is         |
    %|   organized so that it is the variable's elements that are         |
    %|   permuted when the information in the analysis is consumed.       |
    %|   To clarify. Consider som expression 'expr' for which we want to  |
    %|   compute reaching variables. Either we want to test if certain    |
    %|   variables reach the end expression or we want to know which      |
    %|   elements of a variable that is part of the ast reaches the       |
    %|   expression. Say that 'var' reaches the expression and 're' is    |
    %|   our results data. Then 'var(re.reaching_idx)' are the elements   |
    %|   in order that reaches 'expr(re.expr_elem)'. Here 're.expr_elem'  |
    %|   is simply a logical array that picks out the elements of the     |
    %|   expression, while 're.reaching_idx' also does permutation for    |
    %|   the variable.                                                    |
    %|                                                                    |
    %| Properties:                                                        |
    %|   var - The variable the RE analysis concerns.                     |
    %|   enum - Enumeration of all elements of the variable.              |
    %|   reaching - The enumerations that reaches the expression.         |
    %|   expr_elem - The elements of the expression that is reached.      |
    %|                                                                    |
    %| Methods:                                                           |
    %|   enumerate - Enumerate the elements of the variable.              |
    %|   set_enumeration - Set the variables underlying value to the      |
    %|                     enumeration.                                   |
    %|   set_expr_elements_reached - Computes a vector corresponding to   |
    %|                               the elements of the expression the   |
    %|                               enumeration reaches.                 |
    %|   reaching_idx - Returns the index into the variable of the        |
    %|                  reaching elements.                                |
    %|   set_reaching - Computes the indices of the variable the reaches  |
    %|                  the expression. This also preserves ordering of   |
    %|                  the elements.                                     |
    %|____________________________________________________________________|
    properties
        var  % The variable the result concerns
        enum % The enumeration of the elements of the variable
        reaching  % The enumerations that reaches the end expression
        expr_elem % The idxs in the expression that corresponds to the enum
    end
    methods
        function [idx0, idx] = enumerate(obj, idx0)
            %______________________________________________________________
            %|YOP.RE_DATA/ENUMERATE Enumerate all elements.               |
            %|                                                            |
            %| Use:                                                       |
            %|   idx0 = enumerate(obj, idx0)                              |
            %|   idx0 = obj.enumerate(idx0)                               |
            %|   [idx0, idx] = enumerate(obj, idx0)                       |
            %|   [idx0, idx] = obj.enumerate(idx0)                        |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the re_data.                             |
            %|   idx0 - Index to start indexing from.                     |
            %|                                                            |
            %| Return values:                                             |
            %|   idx0 - Incremented start index for next variable.        |
            %|   idx - The enumeration that was used.                     |
            %|____________________________________________________________|
            n = prod(size(obj.var));
            idx = idx0:(idx0+n-1);
            obj.enum = idx;
            obj.set_enumeration();
            idx0 = idx0 + n;
        end
        
        function obj = set_enumeration(obj)
            %______________________________________________________________
            %|YOP.RE_DATA/SET_ENUMERATION Set the value of the underlying |
            %|ast_variable to the enumeration.                            |
            %|                                                            |
            %| Use:                                                       |
            %|   set_enumeration(obj)                                     |
            %|   obj.set_enumeration()                                    |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the re_data.                             |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the re_data.                             |
            %|____________________________________________________________|
            obj.var.m_value = reshape(obj.enum, size(obj.var));
        end
        
        function obj = set_expr_elements_reached(obj, re)
            %______________________________________________________________
            %|YOP.RE_DATA/SET_EXPR_ELEMENTS_REACHED Computes the elements |
            %|in the expression the variable reaches.                     |
            %|                                                            |
            %| Use:                                                       |
            %|   set_expr_elements_reached(obj, re)                       |
            %|   obj.set_expr_elements_reached(re)                        |
            %|                                                            |
            %| Description:                                               |
            %|   Computes a bool array with the same size as the          |
            %|   expression the analysis is made for. A true value        |
            %|   indicate an element the variable reaches.                |
            %|   'expr(re_data(k).expr_elem)' gives the elements of       |
            %|   variable k that reaches the expression.                  |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the re_data.                             |
            %|   re - The result after evaluating the expression.         |
            %|        That is, the reaching elements of all variables.    |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the re_data.                             |
            %|____________________________________________________________|
            obj.set_reaching(re);
            obj.expr_elem = false(size(re));
            for k=1:length(obj.reaching)
                obj.expr_elem(re == obj.reaching(k)) = true;
            end
        end
        
        function idx = reaching_idx(obj)
            %______________________________________________________________
            %|YOP.RE_DATA/REACHING_IDX The indices of the variable that   |
            %|reaches the expression.                                     |
            %|                                                            |
            %| Use:                                                       |
            %|   idx = reaching_idx(obj)                                  |
            %|   idx = obj.reaching_idx()                                 |
            %|                                                            |
            %| Description:                                               |
            %|   The method returns the indices in the order they         |
            %|   reach the expression. For instance                       |
            %|   re_data(k).var(re_data(k).reaching_idx) means that the   |
            %|   the elements of variable k reaches that reaches the      |
            %|   expression are 'reaching_idx' and in that order.         |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the re_data.                             |
            %|                                                            |
            %| Return values:                                             |
            %|   idx - The indices of the variable that reaches the       |
            %|         expression.                                        |
            %|____________________________________________________________|
            idx = obj.reaching - obj.enum(1) + 1;
        end
        
        function obj = set_reaching(obj, re)
            %______________________________________________________________
            %|YOP.RE_DATA/SET_REACHING Pick out the reaching enumeration  |
            %|corresponding to the variable represented by the object.    |
            %|                                                            |
            %| Use:                                                       |
            %|   set_reaching(obj, re)                                    |
            %|   obj.set_reaching(re)                                     |
            %|                                                            |
            %| Parameters:                                                |
            %|   obj - Handle to the re_data.                             |
            %|   re - The enumeration that reaches the expression.        |
            %|                                                            |
            %| Return values:                                             |
            %|   obj - Handle to the re_data.                             |
            %|____________________________________________________________|
            % re - elements that reaches the expression
            obj.reaching = re(re >= obj.enum(1) & re <= obj.enum(end));
        end
        
        function vec = IDs(obj)
            vec = [];
            for k=1:length(obj)
                vec(end+1) = obj(k).var.id;
            end
        end
    end
end