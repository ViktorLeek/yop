function yop_independent(varargin)
% Argument order: t, t0, tf

if nargin >= 1
    name = varargin{1};    
    evalin('caller', ...
        [name ' = yop.ast_independent(''', name, ''');']);
end

if nargin >= 2
    name = varargin{2};    
    evalin('caller', ...
        [name ' = yop.ast_independent_initial(''', name, ''');']);
end

if nargin >=3
    name = varargin{3};    
    evalin('caller', ...
        [name ' = yop.ast_independent_final(''', name, ''');']);
end
    
end