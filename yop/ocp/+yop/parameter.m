function p = parameter(varargin)

ip = inputParser();
ip.FunctionName = "yop.parameter";
ip.addOptional('rows', 1);
ip.addParameter('name', 'p');
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;

if rows == 1
    p = yop.ast_parameter(name);
else
    % This is complex in order to only create one ast_vertcat node
    parameters = cell(rows, 1);
    for k=1:rows
        parameters{k} = yop.ast_parameter([name , '_', num2str(k)]);
    end
    p = vertcat(parameters{:});
end

end