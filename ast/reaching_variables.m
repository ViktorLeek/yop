function res = reaching_variables(expr)
% Get the variables reaching a certain expression
% Enumerates all elements of the variables that make up the expression. It
% den propagates the variables through the expression in order to see which
% reaches the end. Those elements of the final expression that are not
% variables are set to -1. 


% Get the variables that make up the expression
vars = get_variables(expr);

res(length(vars)) = yop.rv_data();

for k=1:length(res)
    res(k).var = vars{k};
end

% Enumerate all elements and set the values of the variables
e0 = 1;
for k=1:length(res)
    e0 = res(k).enumerate(e0); 
end

% Evaluate in order to find reaching variables
reaching_expr = forward_evaluate(expr);

% Elements that are not variables are set to -1
reaching_expr(~isa_variable(expr)) = -1;

% Compute the index in the expression that matches the elements the
% variable takes.
for k=1:length(res)
    res(k).compute_indices(reaching_expr);
end

for k=1:length(res)
    res(k).restore_value();
end

end