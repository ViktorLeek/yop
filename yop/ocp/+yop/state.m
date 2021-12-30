function x = state(varargin)

ip = inputParser();
ip.FunctionName = "yop.state";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'x');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
nx = prod(sz);
name = ip.Results.name;
weight = ones(nx,1) .* ip.Results.scaling(:);
offset = ones(nx,1) .* ip.Results.offset(:);

if isequal(sz, [1, 1])
    x = yop.ast_state(name, weight, offset);
else
    states = cell(nx, 1);
    for k=1:nx
        states{k} = yop.ast_state([name , '_', num2str(k)], weight(k), ...
            offset(k));
    end
    x = reshape(vertcat(states{:}), sz);
end
end