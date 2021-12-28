function x = state(varargin)

ip = inputParser();
ip.FunctionName = "yop.state";
ip.addOptional('rows', 1);
ip.addParameter('name', 'x');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;
weight = ones(rows,1) .* ip.Results.weight(:);
offset = ones(rows,1) .* ip.Results.offset(:);

if rows == 1
    x = yop.ast_state(name, weight, offset);
else
    % This is complex in order to only create one ast_vertcat node
    states = cell(rows, 1);
    for k=1:rows
        states{k} = yop.ast_state([name , '_', num2str(k)], weight(k), ...
            offset(k));
    end
    x = vertcat(states{:});
end
end