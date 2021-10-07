function I = parameterize_integral(int,N,tau,dt,t,x,u,p,tps,ints,ders)
I = 0;
for n=1:N
    yval = [];
    for r=tau
        tt = t(n).evaluate(r);
        xx = x(n).evaluate(r);
        uu = u(n).evaluate(r);
        pp = p;
        val_r = int.fn(tt, xx, uu, pp, tps, ints, ders);
        yval = [yval, val_r(:)];
    end
    lp = yop.lagrange_polynomial(tau, yval).integrate();
    I = I + lp.evaluate(1)*dt;
end
end