function sr = get_subrelation(relation, idx)
% Since indices might be scaled and variables can be scalars it is
% necessary to test if it is possible to take the subindices of the
% relations.

if all(idx==false)
    sr = [];
    return;
end

if isscalar(relation.rhs) %&& ~all(idx==false)
    rhs = relation.rhs;
else
    rhs = relation.rhs(idx);
end

if isscalar(relation.lhs)  %&& ~all(idx==false)
    lhs = relation.lhs;
else
    lhs = relation.lhs(idx);
end

f  = get_constructor(relation);
sr = f(lhs, rhs);

end