function t = independent0(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent0";
ip.addParameter('name', 't0');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

t = yop.ast_independent_initial( ...
    ip.Results.name, ...
    ip.Results.weight, ...
    ip.Results.offset ...
    );

end