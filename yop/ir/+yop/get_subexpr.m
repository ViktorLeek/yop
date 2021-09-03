function se = get_subexpr(expr, idx)

if all(idx==false)
    se = [];
    return;
end

if isscalar(expr)
    se = expr;
else
    se = expr(idx);
end

end