function z = algebraic(varargin)

ip = inputParser();
ip.FunctionName = "yop.algebraic";
ip.addOptional('rows', 1);
ip.addParameter('name', 'z');
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;

if rows == 1
    z = yop.ast_algebraic(name);
else
    % This is complex in order to only create one ast_vertcat node
    algebraics = cell(rows, 1);
    for k=1:rows
        algebraics{k} = yop.ast_algebraic([name , '_', num2str(k)]);
    end
    z = vertcat(algebraics{:});
end
end