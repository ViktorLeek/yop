function t = independent0(name)
switch nargin
    case 0
        t = yop.ast_independent_initial('t0');
    case 1
        t = yop.ast_independent_initial(name);
end
end