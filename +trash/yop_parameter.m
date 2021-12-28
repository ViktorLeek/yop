function yop_parameter(varargin)
for k=1:length(varargin)
    name = varargin{k};
    evalin('caller', ...
        [name ' = yop.ast_parameter(''', name, ''');']);
end
end