clc;
t0 = yop.time0();
tf = yop.timef();
t = yop.time();
x = yop.state(2);
z = yop.algebraic(2);
u = yop.control(2);
p = yop.parameter(2);

sim = yop.simulation();
% sim.add(der(x)==x+z);
% sim.add(x+z==der(x));
% sim.add(x(t0)==[1;2])
% sim.add([1;2]==x(t0));
% sim.add(t0==0);
% sim.add(0==t0);
% sim.add(tf==1);
% sim.add(1==tf);
% sim.add(p == [1;2]);
% sim.add(3 == p);
% sim.add(u - z == 0)
 