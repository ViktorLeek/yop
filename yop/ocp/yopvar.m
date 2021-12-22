function yopvar(varargin)
% t - independent
% t0 - initial value of the independent
% tf - final value of the independent
% x... - anything beginning with x is a state
% z... - anything beginning with z is an algebraic
% u... - anything beginning with u is a control input
% p... - anything beginning with p is a parameter

for k=1:length(varargin)
    name = varargin{k};
    if strcmp(name, 't')
        evalin('caller', ...
        [name ' = yop.independent(''name'', ', '''', name, '''', ');']);
    
    elseif strcmp(name, 't0')
        evalin('caller', ...
        [name ' = yop.independent0(''name'', ', '''', name, '''', ');']);
    
    elseif strcmp(name, 'tf')
        evalin('caller', ...
        [name ' = yop.independentf(''name'', ', '''', name, '''', ');']);
    
    elseif name(1)=='x'
        evalin('caller', ...
        [name ' = yop.state(''name'', ', '''', name, '''', ');']);
    
    elseif name(1)=='z'
        evalin('caller', ...
        [name ' = yop.algebraic(''name'', ', '''', name, '''', ');']);
    
    elseif name(1)=='u'
        evalin('caller', ...
        [name ' = yop.control(''name'', ', '''', name, '''', ');']);
    
    elseif name(1)=='p'
        evalin('caller', ...
        [name ' = yop.parameter(''name'', ', '''', name, '''', ');']);
    
    else 
        error(yop.error.yopvar_not_recognized(name));
    end
end
end