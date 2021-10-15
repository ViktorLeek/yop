function p = parameter(varargin)

ip = inputParser();
ip.FunctionName = "yop.parameter";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'p');
ip.parse(varargin{:});

sz = ip.Results.size;
name = ip.Results.name;

if isscalar(sz)
    sz = [sz, 1];
end

p = yop.ast_parameter(name, sz(1), sz(2));

end