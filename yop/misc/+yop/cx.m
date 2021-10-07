function var = cx(name, rows, cols)
% casadi expression
switch nargin
    case 0
        name = 'v';
        rows = 1;
        cols = 1;
    case 1
        rows = 1;
        cols = 1;
    case 2
        cols = 1;
    otherwise
        % continue
end
% Casadi Expression
if yop.settings.cx_type == yop.settings.MX
    var = casadi.MX.sym(name, rows, cols);
    
elseif yop.settings.cx_type == yop.settings.SX
    var = casadi.SX.sym(name, rows, cols);
    
else
    
end
end