function yop_algebraic(varargin)
for k=1:length(varargin)
    name = varargin{k};
    evalin('caller', ...
        [name ' = yop.ast_algebraic(''', name, ''');']);
end
end