function print(obj)
% Deterimine of argument is empty before calling method.
%
% If this was not in a package, and instead on the path, the method and
% this function could have the same name, and this function would simply
% print '[]' as calls to methods would be resolved differently.

if isempty(obj)
    fprintf('[]');
else
    % Provided that obj is a yop.node, this is not a recursive call, but a 
    % method call.
    print(obj)
end
end