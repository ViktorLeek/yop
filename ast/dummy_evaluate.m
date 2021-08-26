function value = dummy_evaluate(expr)
% Set variables to dummy values in order to propagte values for numerical
% constants.

vars = get_variables(expr);

old_values = cell(size(vars));
for k=1:length(vars)
    old_values{k} = vars{k}.m_value;
    vars{k}.m_value = nan(size(vars{k}));
end

[sort, n] = topological_sort(expr);

for k=1:(n-1)
    forward(sort{k});
end
value = forward(sort{n});

for k=1:length(vars)
    vars{k}.m_value = old_values{k};
end

end