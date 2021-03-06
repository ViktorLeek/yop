function t = independentf(varargin)

ip = inputParser();
ip.FunctionName = "yop.independentf";
ip.addParameter('name', 'tf');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

t = yop.ast_independent_final( ...
    ip.Results.name, ...
    ip.Results.scaling, ...
    ip.Results.offset ...
    );

end