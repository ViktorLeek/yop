function yop_parameter(varargin)
for k=1:length(varargin)
    name = varargin{k};
    % Evaluate in the caller's workspace:
    %   name = yop.state('name', 'name');
    evalin('caller', ...
        [name ' = yop.parameter(''name'', ', '''', name, '''', ');']);
end
end