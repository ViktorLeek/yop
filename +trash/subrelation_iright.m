function sr = subrelation_iright(relation, idx)

assert(numel(idx)==prod(size(relation.right)), '[Yop] Unexpected error');

if all(idx==false)
    sr = [];
    return;
end

rhs = relation.rhs(idx);

if ~isscalar(relation.lhs) && ~isscalar(idx)
    lhs = relation.lhs(idx);
else
    lhs = relation.lhs;
end

f = get_constructor(relation);
sr = f(lhs, rhs);

end