function [t, t0, tf] = independent(name)
switch nargin
    case 0
        t = yop.ast_independent('t');
        t0 = yop.ast_independent_initial('t0');
        tf = yop.ast_independent_final('tf');
        
    case 1
        t = yop.ast_independent(name);
        t0 = yop.ast_independent_initial([name, '0']);
        tf = yop.ast_independent_final([name, 'f']);
end
end