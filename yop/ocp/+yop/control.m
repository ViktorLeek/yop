function x = control(name, rows, cols)
switch nargin
    case 0
        x = yop.ast_control('u', 1, 1);
        
    case 1
        x = yop.ast_control(name, 1, 1);
        
    case 2
        x = yop.ast_control(name, rows, 1);
        
    case 3
        x = yop.ast_control(name, rows, cols);
        
end
end