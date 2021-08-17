function x = state(name, rows, cols)
switch nargin
    case 0
        x = yop.ast_state('x', 1, 1);
        
    case 1
        x = yop.ast_state(name, 1, 1);
        
    case 2
        x = yop.ast_state(name, rows, 1);
        
    case 3
        x = yop.ast_state(name, rows, cols);
        
end
end