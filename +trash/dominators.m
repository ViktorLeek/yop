[tsort, n_elem] = topological_sort(expr);

% Remove nodes that are not yop.nodes.
nodes = tsort(cellfun(@(e) isa(e, 'yop.node'), tsort));
N = length(nodes);

% Reset all predecessors so that only nodes in expr is included
for k=1:N
    reset_pred(nodes{k});
end

% Set predecessors
for k=1:N
    set_pred(nodes{k});
end

% Compute dominators
for k=N:-1:1
    comp_dom(nodes{k}); % Node k is dominated by
end

% Filter out timevarying variables (variables that are parameterized
% differently in the transcription depending on the independent variable)
filt = @(e) ...
    isa(e, 'yop.ast_state') || ...
    isa(e, 'yop.ast_algebraic') || ...
    isa(e, 'yop.ast_control');

vars = nodes(cellfun(filt, nodes));

% If there is a variable that is not dominated by a timepoint or integral,
% then the expression is time_varying
time_varying = true;
tvs = {};
for k=1:length(vars)
    time_invariant_k = false; % time invariant
    for n=1:length(vars{k}.dom)
        dn = vars{k}.dom{n};
        time_invariant_k = time_invariant_k || ...
            isa(dn, 'yop.ast_timepoint') || isa(dn, 'yop.ast_int');
    end
    if ~time_invariant_k
        tvs{end+1} = vars{k};
    end
    time_varying = time_varying && ~time_invariant_k;
end
time_varying



% tf = nodes; % transcription form
% for k=1:length(nodes)
%     nk = nodes{k};
%     if isa(nk, 'yop.ast_int') || isa(nk, 'yop.ast_timepoint')
%         ds = get_ids(nk.dominates);
%         ds = ds(2:end); % every node dominates itself.
%         isec = intersect(get_ids(tf), ds);
%     end
% end

% Compute immediate dominators
% for k=1:N
%     if isa(nodes{k}, 'yop.node')
%         nodes{k}.idom = sdom(nodes{k});
%     end
% end
% for n=(N-1):-1:1 % not root
%     idom = nodes{n}.idom;
%     tmp = idom;
%     idom_n = 1:length(tmp);
%     for s=idom_n
%         for t=idom_n(idom_n~=s)
%             if ~isempty( find(yop.get_ids(tmp{s}.idom)==tmp{t}.id, 1) )
%                 %t in tsort{s}.idom
%                 idom{t} = [];
%             end
%         end
%     end
%     nodes{n}.idom = idom(~cellfun('isempty', idom));
% end

% Traverse and set immediately dominates property in order to obtain a dom
% tree

% Traverse tree in order to evaluate or new tsort