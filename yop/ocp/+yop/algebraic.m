function z = algebraic(varargin)

ip = inputParser();
ip.FunctionName = "yop.algebraic";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'z');
ip.parse(varargin{:});

sz = ip.Results.size;
name = ip.Results.name;

if isscalar(sz)
    sz = [sz, 1];
end

z = yop.ast_algebraic(name, sz(1), sz(2));

end