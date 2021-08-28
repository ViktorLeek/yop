function vars = get_variables(expressions)
vars = {};

if ~isa(expressions, 'cell')
    expressions = {expressions};
end

% First iteration is hoisted out in order to later warm start the search
% for every expression. Otherwise it might even lead to erroneous 
% results as doublettes are likely to appear, these could be found 
% afterwards, but it is faster to warm start the search after the 
% first iteration.

[tsort, n_elem, visited] = topological_sort(expressions{1});
for n=1:n_elem
    if isa(tsort{n}, 'yop.ast_variable')
        vars = {vars{:}, tsort{n}};
    end
end

for k=2:length(expressions)
    
    [tsort, n_elem, visited] = topological_sort(expressions{k}, visited);
    
    for n=1:n_elem
        if isa(tsort{n}, 'yop.ast_variable')
            vars = {vars{:}, tsort{n}};
        end
    end
end
end