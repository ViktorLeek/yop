function z = algebraic(varargin)

ip = inputParser();
ip.FunctionName = "yop.algebraic";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'z');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
nz = prod(sz);
name = ip.Results.name;
weight = ones(nz,1) .* ip.Results.weight(:);
offset = ones(nz,1) .* ip.Results.offset(:);

if isequal(sz, [1, 1])
    z = yop.ast_algebraic(name, weight, offset);
else
    algs = cell(nz, 1);
    for k=1:nz
        algs{k} = yop.ast_algebraic([name , '_', num2str(k)], weight(k), ...
            offset(k));
    end
    z = reshape(vertcat(algs{:}), sz);
end
end