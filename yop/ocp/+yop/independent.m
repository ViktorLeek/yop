function t = independent(name)
switch nargin
    case 0
        t = yop.ast_independent('t');
    case 1
        t = yop.ast_independent(name);
end
end