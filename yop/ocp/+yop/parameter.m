function p = parameter(varargin)

ip = inputParser();
ip.FunctionName = "yop.parameter";
ip.addOptional('rows', 1);
ip.addParameter('name', 'p');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;
weight = ones(rows,1) .* ip.Results.weight(:);
offset = ones(rows,1) .* ip.Results.offset(:);

if rows == 1
    p = yop.ast_parameter(name, weight, offset);
else
    % This is complex in order to only create one ast_vertcat node
    parameters = cell(rows, 1);
    for k=1:rows
        parameters{k} = yop.ast_parameter([name , '_', num2str(k)], ...
            weight(k), offset(k));
    end
    p = vertcat(parameters{:});
end

end