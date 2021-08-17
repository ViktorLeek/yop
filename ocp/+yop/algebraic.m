function z = algebraic(name, rows, cols)
switch nargin
    case 0
        z = yop.ast_algebraic('z', 1, 1);
        
    case 1
        z = yop.ast_algebraic(name, 1, 1);
        
    case 2
        z = yop.ast_algebraic(name, rows, 1);
        
    case 3
        z = yop.ast_algebraic(name, rows, cols);
        
end
end