function p = parameter(varargin)

ip = inputParser();
ip.FunctionName = "yop.parameter";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'p');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
np = prod(sz);
name = ip.Results.name;
weight = ones(np,1) .* ip.Results.weight(:);
offset = ones(np,1) .* ip.Results.offset(:);

if isequal(sz, [1, 1])
    p = yop.ast_parameter(name, weight, offset);
else
    params = cell(np, 1);
    for k=1:np
        params{k} = yop.ast_parameter([name , '_', num2str(k)], weight(k), ...
            offset(k));
    end
    p = reshape(vertcat(params{:}), sz);
end
end