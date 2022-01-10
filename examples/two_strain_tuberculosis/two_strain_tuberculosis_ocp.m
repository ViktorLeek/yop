% yops times: t t0 tf
% yops state: x size: [6,1]
% yops ctrl:  u size: [2,1]
% 
% N_pop = 30000;
% x0   = [76; 1; 36; 2; 4; 1]*Npop/120;
% 
% R = diag([25; 250]);
% ocp = yop.ocp('Two-Strain Tuberculosis');
% ocp.min( int(x(4) + x(5) + u'*R*u) )
% ocp.st( tf == 5 );
% ocp.st( der(x) == tuberculosis(x, u) );
% ocp.st(  x(t0) == x0 );
% ocp.st( 0.00 <= x <= 30e3 );
% ocp.st( 0.05 <= u <= 0.95 );
% ocp.st( sum(x) == N_pop );
% 
% timeGuess = [t0; tf];
% SGuess
% TGuess
% L1Guess
% L2Guess
% I1Guess
% I2Guess
% u1Guess
% u2Guess
% guess.phase.time    = [timeGuess];
% guess.phase.state   = [SGuess, TGuess, L1Guess, L2Guess, I1Guess, I2Guess];
% guess.phase.control = [u1Guess, u2Guess];
% guess.phase.integral = 6000;

ig = yop.guess(t, tval, x, xval, z, zval, u, uval, p, pval);