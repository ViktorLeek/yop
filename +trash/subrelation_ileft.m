function sr = subrelation_ileft(relation, idx)

assert(numel(idx)==prod(size(relation.left)), '[Yop] Unexpected error');

if all(idx==false)
    sr = [];
    return;
end

lhs = relation.lhs(idx);

if ~isscalar(relation.rhs) && ~isscalar(idx)
    rhs = relation.rhs(idx);
else
    rhs = relation.rhs;
end

f = get_constructor(relation);
sr = f(lhs, rhs);

end