function [re, nr, reaching_elements] = reaching_elems(expr, vars)
% re = reaching_elems(expr)
% Get the elements of a variable that reaches a certain expression.
%   expr - The expression the elements reach
%   re   - A vector of yop.re_data with one element per variable that
%          reaches the expression. Every element has four fields:
%          1) .var      - A variable found in the expression tree.
%          2) .enum     - An enumeration of all elements of that variable.
%                         If there are more than one variable, indexation
%                         starts at one for the first variable, and then
%                         increases by one for every element.
%          3) .reaching - The elements that reaches the expression.
%          4) .expr_elem - A logical array the same size as the expression. 
%                         True values are those elements that reaches the
%                         final expression. That is, expr(res(k).index)
%                         gives the elements of variable k that reaches the
%                         expression.
%          5) .reaching_idx  - The indices of the variable that reaches the 
%                         expression. 'var_k(res(k).reaching_idx)' returns the
%                         elements of the variable that reaches the
%                         expression. The benefit of this is if the
%                         enumeration is 5:8, and [8,5] reaches the
%                         expression (in that order). Then 
%                         var(reaching_idx) == expr(expr_elem).


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
reaching_elements = propagate_value(expr);

% Elements that are not variables are set to 0
reaching_elements(~isa_variable(expr)) = 0;

% Compute the index in the expression that matches the elements the
% variable takes.
for k=1:length(re)
    re(k).set_expr_elements_reached(reaching_elements);
end

% Remove variables that do not reach the final expression.
reaching = arrayfun(@(v) ~isempty(v.reaching), re);

nr = re(~reaching);
re = re( reaching);

end