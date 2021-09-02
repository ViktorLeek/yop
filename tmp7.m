%% Expression
[t, t0, tf] = yop.time();
p = yop.parameter();
x1 = yop.state();
v = [x1, x1+p];
v(1) = x1(t==1);
expr = v;
is_transcription_invariant(expr)

