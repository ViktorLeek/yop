function vars = get_variables(expressions)
vars = {};

% First iteration is hoisted out in order to later warm start the search
% for every expression. Otherwise it might even lead to erroneous 
% results as doublettes are likely to appear, these could be found 
% afterwards, but it is faster to warm start the search after the 
% first iteration.

[tsort, visited, n_elem] = topological_sort(expressions{1});
for n=1:n_elem
    if isa(tsort{n}, 'yop.ast_variable')
        % For small cell arrays, measurements indicate that creating a new
        % array seems to be faster, as opposed to Matlab's suggestions.
        vars = {vars{:}, tsort{n}};
        % vars = [vars(:)', tsort(n)'];
    end
end

for k=2:length(expressions)
    % By preallocation the sort is sped up considerably
    tmp = cell(yop.constants().topsort_preallocation_size, 1);
    
    [tsort, visited, n_elem] = topological_sort( ...
        expressions{k}, ... Warm start procedure below:
        tmp,            ...  - Nodes sorted so far
        visited,        ...  - Warm start, only visit new nodes
        0               ...  - Number of nodes sorted
        );
    
    for n=1:n_elem
        if isa(tsort{n}, 'yop.ast_variable')
            vars = {vars{:}, tsort{n}};
            % vars = [vars(:)', tsort(n)'];
        end
    end
end
end