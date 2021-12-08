function x = state(varargin)

ip = inputParser();
ip.FunctionName = "yop.state";
ip.addOptional('rows', 1);
ip.addParameter('name', 'x');
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;

if rows == 1
    x = yop.ast_state(name);
else
    % This is complex in order to only create one ast_vertcat node
    states = cell(rows, 1);
    for k=1:rows
        states{k} = yop.ast_state([name , '_', num2str(k)]);
    end
    x = vertcat(states{:});
end
end