% PROPT 44
toms t t_f

% Parameters
aalpha = 0.01227; bbeta = 0.145e-3;

c  = 2060;    g0 = 9.81;
r0 = 6.371e6; r02=r0*r0;
m0 = 215; 
mf = 68;
Fm = 9.525515;

n=60;

    
p = tomPhase('p', t, 0, t_f, n);
setPhase(p);
tomStates h v m
tomControls F

% Initial guess

x0 = { ...
    t_f == tfopt, ...
    icollocate({v == vopt; h == hopt; m == mopt})
    collocate(F == Fopt) ...
    };

% Box constraints
cbox = { ...
    100 <= t_f <= 300,
    icollocate({0 <= v; 0 <= h; mf <= m <= m0; 0 <= F <= Fm}) ...
    };

% Boundary constraints
cbnd = {initial({v == 0; h == 0; m == m0})
    final({v==0; m == mf})};

D = aalpha*v.^2.*exp(-bbeta*h);
g = g0; % or g0*r02./(r0+h).^2;

% ODEs and path constraints
ceq = collocate({dot(h) == v, m*dot(v) == F*c-D-g*m, dot(m) == -F});

% Objective
objective = -1e-4*final(h);

options = struct;
options.name = 'Goddard Rocket 1';
options.Prob.SOL.optPar(30) = 30000;
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);

% Optimal v and more to use as starting guess
vopt = subs(v, solution);
hopt = subs(h, solution);
mopt = subs(m, solution);
Fopt = subs(F, solution);
tfopt = subs(t_f, solution);
    

t = subs(collocate(t),solution);
v = subs(collocate(vopt),solution);
h = subs(collocate(hopt),solution);
m = subs(collocate(mopt),solution);
F = subs(collocate(Fopt),solution);

figure(1)
subplot(311); hold on
plot(t, v)
xlabel('Time'); ylabel('Velocity')

subplot(312); hold on
plot(t, h)
xlabel('Time'); ylabel('Height')

subplot(313); hold on
plot(t, m)
xlabel('Time'); ylabel('Mass')

figure(2); hold on
stairs(t, F)
xlabel('Time'); ylabel('F (Control)')