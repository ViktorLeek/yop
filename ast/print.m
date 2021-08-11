function print(obj)
% print(obj) - Part of Yop
% Gives polymorphic behaviour to built in classes from an AST perspective.
% The function needs to be on the path and outside the package in order to
% be polymorphic.

if isempty(obj)
    fprintf('[]');
    
elseif isnumeric(obj)
    % By making several calls to fprintf, potential dimensionality problems
    % are avoided when concatenating strings.
    fprintf('[');
    fprintf(num2str(obj)); 
    fprintf(']\n');
    
else
    % Unknown type, try to print.
    fprintf(obj); 
    fprintf('\n');
end
end