function t = independent(varargin)

ip = inputParser();
ip.FunctionName = "yop.independent";
ip.addParameter('name', 't');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

t = yop.ast_independent( ...
    ip.Results.name, ...
    ip.Results.scaling, ...
    ip.Results.offset ...
    );

end