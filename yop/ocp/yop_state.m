function yop_state(varargin)
for k=1:length(varargin)
    name = varargin{k};
    % Evaluate in the caller's workspace:
    %   name = yop.state('name', 'name');
    evalin('caller', ...
        [name ' = yop.state(''name'', ', '''', name, '''', ');']);
end
end