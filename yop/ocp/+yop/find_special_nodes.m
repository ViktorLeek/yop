function [vars, tps, ints, ders, sorted] = find_special_nodes(expr)

vars = {};
sorted = yop.ocp_expr.empty(1,0); % Topological order
tps   = yop.ocp_expr.empty(1,0);
ints  = yop.ocp_expr.empty(1,0);
ders  = yop.ocp_expr.empty(1,0);

% Sort all nodes
[tsort, n_elem] = topological_sort(expr);

% Find special nodes, maintin topological order in sorted.
for n=1:n_elem
    if isa(tsort{n}, 'yop.ast_variable')
        vars{end+1} = tsort{n};
        
    elseif isa(tsort{n}, 'yop.ast_timepoint')
        sn = yop.ocp_expr(tsort{n}, yop.ocp_expr.tp);
        tps(end+1) = sn;
        sorted(end+1) = sn;
        
    elseif isa(tsort{n}, 'yop.ast_int')
        sn = yop.ocp_expr(tsort{n}, yop.ocp_expr.int);
        ints(end+1) = sn;
        sorted(end+1) = sn;
        
    elseif isa(tsort{n}, 'yop.ast_der')
        sn = yop.ocp_expr(tsort{n}, yop.ocp_expr.der);
        ders(end+1) = sn;
        sorted(end+1) = sn;
        
    end
end
end