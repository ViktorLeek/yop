function t = independent(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent";
ip.addParameter('name', 't');
ip.parse(varargin{:});

t = yop.ast_independent(ip.Results.name);

end