function srf = to_srf(constraints)
% to_srl - To single relation form

% Find all 'ast_relation' nodes in the constraints
relations = {};
for k=1:length(constraints)
    relations = [relations(:)', get_relations(constraints{k})];
end

% Convert the relations to srf
srf = {};
for k=1:length(relations)
    rk = relations{k};
    switch class(rk)
        case 'yop.ast_lt'
            srf = {srf{:}, yop.ast_lt(rmost(rk.lhs), lmost(rk.rhs))};
            
        case 'yop.ast_gt'
            srf = {srf{:}, yop.ast_gt(rmost(rk.lhs), lmost(rk.rhs))};
            
        case 'yop.ast_le'
            srf = {srf{:}, yop.ast_le(rmost(rk.lhs), lmost(rk.rhs))};
            
        case 'yop.ast_ge'
            srf = {srf{:}, yop.ast_ge(rmost(rk.lhs), lmost(rk.rhs))};
            
        case 'yop.ast_eq'
            srf = {srf{:}, yop.ast_eq(rmost(rk.lhs), lmost(rk.rhs))};
            
        case 'yop.ast_ne'            
            srf = {srf{:}, yop.ast_ne(rmost(rk.lhs), lmost(rk.rhs))};
            
        otherwise
            error('[yop] Error: unknown relation')
    end
end

end