function t = independent0(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent0";
ip.addParameter('name', 't0');
ip.parse(varargin{:});

name = ip.Results.name;
t = yop.ast_independent_initial(name);

end