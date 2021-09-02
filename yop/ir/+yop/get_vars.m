function vars = get_vars(expressions_cell)
% get variables of expressions in cell array.
% A bit of a stupid name, because there is a yop.node method with the same
% name that operates on individual nodes. Because of this, this function is
% located withing the yop namespace, so to call it use yop.get_vars.

vars = {};

if ~isa(expressions_cell, 'cell')
    expressions_cell = {expressions_cell};
end

% First iteration is hoisted out in order to later warm start the search
% for every expression. Otherwise it might even lead to erroneous 
% results as doublettes are likely to appear, these could be found 
% afterwards, but it is faster to warm start the search after the 
% first iteration.

[tsort, n_elem, visited] = topological_sort(expressions_cell{1});
for n=1:n_elem
    if isa(tsort{n}, 'yop.ast_variable')
        vars = {vars{:}, tsort{n}};
    end
end

for k=2:length(expressions_cell)
    
    [tsort, n_elem, visited] = topological_sort(expressions_cell{k}, visited);
    
    for n=1:n_elem
        if isa(tsort{n}, 'yop.ast_variable')
            vars = {vars{:}, tsort{n}};
        end
    end
end
end