function draw(obj)
% ast - print the abstract syntax tree of a yop expression.

% This file gives polymorphic behaviour to built in classes from an AST 
% perspective. % The function needs to be on the path and outside the 
% package in order to be polymorphic.

if isempty(obj)
    fprintf('[]\n');
    
elseif isnumeric(obj)
    % By making several calls to fprintf, potential dimensionality problems
    % are avoided when concatenating strings.
    fprintf('[');
    fprintf(num2str(obj)); 
    fprintf(']\n');
    
elseif isstring(obj)
    fprintf('\"');
    fprintf(obj);
    fprintf('\"\n');
    
elseif ischar(obj)
    fprintf([char(39), obj, char(39), '\n']);
    
else
    % Unknown type, try to print.
    fprintf(obj); 
    fprintf('\n');
end

end