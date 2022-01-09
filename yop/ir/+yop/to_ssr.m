function ssr = to_ssr(constraint)
% To scalar single relation form
relations = get_relations(constraint);
K = length(relations);
ssr = {};
for k=1:K
    rk = relations{k};
    rel = get_constructor(rk);
    sr = rel(rmost(rk.m_lhs), lmost(rk.m_rhs));
    numel_lhs = prod(size(sr.m_lhs));
    numel_rhs = prod(size(sr.m_rhs));
    if numel_lhs > 1 && numel_rhs == 1
        for n=1:numel_lhs
            rhs = get_subexpr(sr.m_rhs, n); % size(@fn) == [1,1]
            ssr{end+1} = rel(sr.m_lhs(n), rhs);
        end
    elseif numel_lhs==1 && numel_rhs > 1
        for n=1:numel_rhs
            lhs = get_subexpr(sr.m_lhs, n); % size(@fn) == [1,1]
            ssr{end+1} = rel(lhs, sr.m_rhs(n));
        end
    elseif numel_lhs == numel_rhs
        for n=1:numel_lhs
            lhs = get_subexpr(sr.m_lhs, n); % size(@fn) == [1,1]
            rhs = get_subexpr(sr.m_rhs, n); % size(@fn) == [1,1]
            ssr{end+1} = rel(lhs, rhs);
        end
    else
        error(yop.error.incompatible_constraint_size());
    end
end
end

function sexpr = get_subexpr(expr, n)
if isscalar(expr) && isnumeric(expr)
    sexpr = expr;
elseif isa(expr, 'function_handle')
    sexpr = @(t) yop.get_elem(expr(t), n);
else
    sexpr = expr(n);
end
end