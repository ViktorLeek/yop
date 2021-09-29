function val = parameterize_timepoint(tp, t0, tf, t, x, u, p, tps, ints, ders)
tt = t.value(tp.timepoint, t0, tf);
xx = x.value(tp.timepoint, t0, tf);
uu = u.value(tp.timepoint, t0, tf); 
val = tp.fn(tt, xx, uu, p, tps, ints, ders);
end