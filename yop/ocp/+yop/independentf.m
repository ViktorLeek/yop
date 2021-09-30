function t = independentf(name)
switch nargin
    case 0
        t = yop.ast_independent_final('tf');
    case 1
        t = yop.ast_independent_final(name);
end
end