function u = control(name, rows, cols)
% Ändra till varargout. Antal variabler beror på pwx.
switch nargin
    case 0
        u = yop.ast_control('u', 1, 1);
        
    case 1
        u = yop.ast_control(name, 1, 1);
        
    case 2
        u = yop.ast_control(name, rows, 1);
        
    case 3
        u = yop.ast_control(name, rows, cols);
        
    case 4
        u = yop.ast_control(name, rows, cols);
        
end
end