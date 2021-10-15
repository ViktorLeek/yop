function x = state(varargin)

ip = inputParser();
ip.FunctionName = "yop.state";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'x');
ip.parse(varargin{:});

sz = ip.Results.size;
name = ip.Results.name;

if isscalar(sz)
    sz = [sz, 1];
end

x = yop.ast_state(name, sz(1), sz(2));

end