function p = parameter(varargin)

ip = inputParser();
ip.FunctionName = "yop.parameter";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'p');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
np = prod(sz);
name = ip.Results.name;
scaling = ones(np,1) .* ip.Results.scaling(:);
offset = ones(np,1) .* ip.Results.offset(:);

if isequal(sz, [1, 1])
    p = yop.ast_parameter(name, scaling, offset);
else
    params = cell(np, 1);
    for k=1:np
        params{k} = yop.ast_parameter([name , '_', num2str(k)], scaling(k), ...
            offset(k));
    end
    p = reshape(vertcat(params{:}), sz);
end
end