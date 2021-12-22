function yop_algebraic(varargin)
for k=1:length(varargin)
    name = varargin{k};
    % Evaluate in the caller's workspace:
    %   name = yop.algebraic('name', 'name');
    evalin('caller', ...
        [name ' = yop.algebraic(''name'', ', '''', name, '''', ');']);
end
end