function yop_alg(varargin)
for k=1:length(varargin)
    name = varargin{k};
    evalin('caller', ...
        [name ' = yop.ast_algebraic(''', name, ''');']);
end
end