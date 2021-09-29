function disc = parameterize_expression(expr,N,tau,t,x,u,p,tps,ints,ders)
if is_transcription_invariant(expr)
    tt = t(1).evaluate(0);
    xx = x(1).evaluate(0);
    uu = u(1).evaluate(0);
    pp = p;
    % dd = der(1).evaluate(0);
    disc = expr.fn(tt, xx, uu, pp, tps, ints, ders);
else
    disc = [];
    for n=1:N
        tt = t(n).evaluate(0);
        xx = x(n).evaluate(0);
        uu = u(n).evaluate(0);
        pp = p;
        % dd = der(n).evaluate(0);
        disc = [disc; expr.fn(tt, xx, uu, pp, tps, ints, ders)];
        if expr.is_hard
            for r=tau(2:end)
                tt = t(n).evaluate(r);
                xx = x(n).evaluate(r);
                uu = u(n).evaluate(r);
                disc = [disc; expr.fn(tt, xx, uu, pp, tps, ints, ders)];
            end
        end
    end
    tt = t(N+1).evaluate(0);
    xx = x(N+1).evaluate(0);
    uu = u(N+1).evaluate(0);
    pp = p;
    disc = [disc; expr.fn(tt, xx, uu, pp, tps, ints, ders)];
end
end