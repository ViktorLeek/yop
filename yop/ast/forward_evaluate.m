function value = forward_evaluate(expr)

[sort, n] = topological_sort(expr);

for k=1:(n-1)
    forward(sort{k});
end
value = forward(sort{n});

end