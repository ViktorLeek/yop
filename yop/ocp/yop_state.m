function yop_state(varargin)
for k=1:length(varargin)
    name = varargin{k};
    evalin('caller', ...
        [name ' = yop.ast_state(''', name, ''');']);
end
end