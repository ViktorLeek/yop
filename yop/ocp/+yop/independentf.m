function t = independentf(varargin)

ip = inputParser();
ip.FunctionName = "yop.independentf";
ip.addParameter('name', 'tf');
ip.parse(varargin{:});

name = ip.Results.name;
t = yop.ast_independent_final(name);

end