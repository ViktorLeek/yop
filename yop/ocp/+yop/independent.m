function t = independent(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent";
ip.addParameter('name', 't');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

t = yop.ast_independent( ...
    ip.Results.name, ...
    ip.Results.weight, ...
    ip.Results.offset ...
    );

end