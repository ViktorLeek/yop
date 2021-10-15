function t = independent(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent";
ip.addParameter('name', 't');
ip.parse(varargin{:});

name = ip.Results.name;
t = yop.ast_independent(name);

end