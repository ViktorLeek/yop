function [re, nr, reaching_enumeration] = reaching_elems(expr, vars)
%__________________________________________________________________________
%|YOP.REACHING_ELEMS Reching elements analysis. Test which elements of    |
%|the variables the is part of the expression, or could be part of the    |
%|expression, that reaches the final expression. An element is reaching   |
%|the expression if it is possible to deterime that the entrie in the     |
%|expression is a variable, which essentially means that it is only       |
%|permuted.                                                               |
%|                                                                        |
%| Use:                                                                   |
%|   re = yop.reaching_elems(expr)                                        |
%|   [re, nr] = yop.reaching_elems(expr, vars)                            |
%|   [re, nr, reaching_enumeration] = yop.reaching_elems(expr, vars)      |
%|                                                                        |
%| Parameters:                                                            |
%|   expr - Ast node for the expression of interest.                      |
%|   vars - (Optional) Variable for which it is of interest to see if     |
%|          and which elements reaches the expression. If this parameter  |
%|          is left out Yop will detect the variables that are part of    |
%|          the expression and do the analysis on them.                   |
%|                                                                        |
%| Return values:                                                         |
%|   re - A vector of yop.re_data, one element for every variable         |
%|        that is reaching the expression.                                |
%|   nr - A vector of yop.re_data, one element for every variable         |
%|        that is not reaching the expression.                            |
%|   reaching_enumeration - The enumeration that is reaching the          |
%|                          expression. Zeroes for those that are not     |
%|                          variables.                                    |
%|________________________________________________________________________|

% Get the variables that make up the expression
if nargin == 1
    % This way it is possible to include, manually, the variables that you
    % want to test if they reach the expression. For instance to answer the
    % question of how many of the state variables the reached the ode rhs
    % expression. Note that the manual vars class should include the
    % get_vars class for proper analysis.
    vars = yop.get_vars(expr);
end

re(length(vars)) = yop.re_data();

for k=1:length(re)
    re(k).var = vars{k};
end

% Enumerate all elements and set the values of the variables
e0 = 1;
for k=1:length(re)
    e0 = re(k).enumerate(e0); 
end

% Evaluate in order to find reaching variables
reaching_enumeration = propagate_value(expr);

% Elements that are not variables are set to 0
reaching_enumeration(~isa_variable(expr)) = 0;

% Compute the index in the expression that matches the elements the
% variable takes.
for k=1:length(re)
    re(k).set_expr_elements_reached(reaching_enumeration);
end

% Remove variables that do not reach the final expression.
reaching = arrayfun(@(v) ~isempty(v.reaching), re);

nr = re(~reaching);
re = re( reaching);

end