function ssr = to_ssr(constraint)
% To scalar single relation form
relations = get_relations(constraint);
K = length(relations);
ssr = {};
for k=1:K
    rk = relations{k};
    rel = get_constructor(rk);
    sr = rel(rmost(rk.lhs), lmost(rk.rhs));
    numel_lhs = prod(size(sr.lhs));
    numel_rhs = prod(size(sr.rhs));
    if numel_lhs > 1 && numel_rhs == 1
        for n=1:numel_lhs
            ssr{end+1} = rel(sr.lhs(n), sr.rhs);
        end
    elseif numel_lhs==1 && numel_rhs > 1
        for n=1:numel_rhs
            ssr{end+1} = rel(sr.lhs, sr.rhs(n));
        end
    elseif numel_lhs == numel_rhs
        for n=1:numel_lhs
            ssr{end+1} = rel(sr.lhs(n), sr.rhs(n));
        end
    else
        error(yop.error.incompatible_constraint_size());
    end
end
end