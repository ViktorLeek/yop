function p = parameter(name, rows, cols)
switch nargin
    case 0
        p = yop.ast_parameter('p', 1, 1);
        
    case 1
        p = yop.ast_parameter(name, 1, 1);
        
    case 2
        p = yop.ast_parameter(name, rows, 1);
        
    case 3
        p = yop.ast_parameter(name, rows, cols);
        
end
end