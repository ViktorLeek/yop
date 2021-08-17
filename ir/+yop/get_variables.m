function vars = get_variables(expressions)
vars = {};
for k=1:length(expressions)
    tsort = toplogical_sort(expressions{k});
    for n=1:length(tsort)
        if isa(tsort{n}, 'yop.ast_variable')
            vars = [vars(:)', tsort(n)'];
        end
    end
end
end