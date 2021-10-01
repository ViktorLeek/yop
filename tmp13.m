% y = yop.algebraic()
% y = yop.external();
y = yop.output();
u = yop.control();
y_ref = yop.control();
x = yop.state();
I = yop.state();
dy = yop.state();

     der(x) == f(t, x, z, u, p);
          y == g(t, x, z, u, p);
[der(I); u] == pid(I, y, y_ref, der(y));
      y_ref == step(0, 1).at(t==1);
      
      p_pid = [Kp; Ki; 0];
      
          t0 == 0;
          tf == T;
      der(x) == f(t, x, u, p);
          y  == o(t, x, u, p);
         yr  == step(0, 1).at(t==1);
          e  == yr-y;
          u  == PID.u(I, e, p_pid);
      der(I) == PID.I(I, e)
      
      
       
      