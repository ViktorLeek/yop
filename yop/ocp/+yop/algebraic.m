function z = algebraic(varargin)

ip = inputParser();
ip.FunctionName = "yop.algebraic";
ip.addOptional('rows', 1);
ip.addParameter('name', 'z');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;
weight = ones(rows,1) .* ip.Results.weight(:);
offset = ones(rows,1) .* ip.Results.offset(:);

if rows == 1
    z = yop.ast_algebraic(name, weight, offset);
else
    algebraics = cell(rows, 1);
    for k=1:rows
        algebraics{k} = yop.ast_algebraic([name , '_', num2str(k)], ...
            weight(k), offset(k));
    end
    z = vertcat(algebraics{:}); % single ast node
end
end