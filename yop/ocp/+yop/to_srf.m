function srf = to_srf(constraints_cell, srf)
% to_srl - To single relation form
%   From the current node, finds all relations, and creates a
%   node with a single relation, for all relations withing this
%   node.

if nargin == 1
    srf = yop.srf_data();
end

for n=1:length(constraints_cell)
    % Find all 'ast_relation' nodes in the constraints
    relations = get_relations(constraints_cell{n});
    
    % Convert the relations to srf
    for k=1:length(relations)
        rk = relations{k};
        switch class(rk)
            case 'yop.ast_lt'
                srf.add_lt(yop.ast_lt(rmost(rk.lhs), lmost(rk.rhs)));
                
            case 'yop.ast_gt'
                srf.add_gt(yop.ast_gt(rmost(rk.lhs), lmost(rk.rhs)));
                
            case 'yop.ast_le'
                srf.add_le(yop.ast_le(rmost(rk.lhs), lmost(rk.rhs)));
                
            case 'yop.ast_ge'
                srf.add_ge(yop.ast_ge(rmost(rk.lhs), lmost(rk.rhs)));
                
            case 'yop.ast_eq'
                srf.add_eq(yop.ast_eq(rmost(rk.lhs), lmost(rk.rhs)));
                
            case 'yop.ast_ne'
                srf.add_ne(yop.ast_ne(rmost(rk.lhs), lmost(rk.rhs)));
                
            otherwise
                error('[yop] Error: unknown relation')
        end
    end
end

end