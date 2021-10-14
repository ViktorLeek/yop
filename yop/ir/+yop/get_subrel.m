function sr = get_subrel(relation, idx)
%GET_SUBREL - Get subrelation from index.
%  sr = yop.get_subrelation(relation, idx)
%
%  Description:
%    Get the subrelations of a relation based on an index. Since the lhs
%    and rhs can have different size it is necessary to test if one side is
%    a scalar and in that case use the entire side in the subrelation.
%
%  Parameters:
%    relation - A relation that meets the specifiction
%               isa(relation, 'yop.ast_relation')
%    idx      - Indices to extract
%               Type: logical
%
%  Example:
%    sr = yop.get_subrelation(relation, isa_variable(relation.lhs))

assert(islogical(idx), '[Yop] Unexpected error.')

if all(idx==false)
    sr = [];
    return;
end

if isscalar(idx)
    lhs = relation.lhs;
    rhs = relation.rhs;
    
elseif isscalar(relation.lhs) && ~isscalar(relation.rhs)
    lhs = relation.lhs;
    rhs = relation.rhs(idx);
    
elseif ~isscalar(relation.lhs) && isscalar(relation.rhs)
    lhs = relation.lhs(idx);
    rhs = relation.rhs;
    
elseif ~isscalar(relation.lhs) && ~isscalar(relation.rhs)
    lhs = relation.lhs(idx);
    rhs = relation.rhs(idx);
    
end

if isa(lhs, 'function_handle')
    lhs = @(t) yop.get_subexpr(lhs(t), idx);
end

if isa(rhs, 'function_handle')
    rhs = @(t) yop.get_subexpr(rhs(t), idx);
end

f = get_constructor(relation);
sr = f(lhs, rhs);
end