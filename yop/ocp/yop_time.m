function yop_time(varargin)
% Argument order: t, t0, tf

if nargin >= 1
    name = varargin{1};    
    evalin('caller', ...
        [name ' = yop.independent(''name'', ', '''', name, '''', ');']);
end

if nargin >= 2
    name = varargin{2};    
    evalin('caller', ...
        [name ' = yop.independent0(''name'', ', '''', name, '''', ');']);
end

if nargin >=3
    name = varargin{3};    
    evalin('caller', ...
        [name ' = yop.independentf(''name'', ', '''', name, '''', ');']);
end
    
end