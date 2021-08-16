function srf = to_srf(constraints)
% to_srl - To single relation form

% Get the relations on the constraints
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
            srf = [srf(:)', { yop.srf_lt(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        case 'yop.ast_gt'
            srf = [srf(:)', { yop.srf_gt(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        case 'yop.ast_le'
            srf = [srf(:)', { yop.srf_le(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        case 'yop.ast_ge'
            srf = [srf(:)', { yop.srf_ge(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        case 'yop.ast_eq'
            srf = [srf(:)', { yop.srf_eq(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        case 'yop.ast_ne'
            srf = [srf(:)', { yop.srf_ne(rmost(rk.lhs), lmost(rk.rhs)) }];
            
        otherwise
            error('[yop] Error: unknown relation')
    end
end

end