tmp = 1:3;

fn = @(t) get_subexpr(tmp, [false, true, false]);
fn(2)


function value = get_subexpr(expr, sub)
value = expr(sub);
end