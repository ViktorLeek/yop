function varargout = control(varargin)

ip = inputParser();
ip.FunctionName = "yop.control";
ip.addOptional('rows', 1);
ip.addParameter('name', 'u');
ip.addParameter('weight', 1);
ip.addParameter('offset', 0);
ip.addParameter('deg', 0);
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;
weight = ones(rows,1) .* ip.Results.weight(:);
offset = ones(rows,1) .* ip.Results.offset(:);
deg = ones(rows,1) .* ip.Results.deg;

if rows <= 0
    error(yop.error.variable_too_few_dimensions());
end

if any(deg < 0)
    error(yop.error.invalid_control_degree());
end

if rows == 1
    u = yop.ast_control(name, weight, offset, deg);
    varargout = {u};
    du = u.der;
    while ~isempty(du)
        varargout{end+1} = du;
        du = du.der;
    end
else
   controls = {}; 
    for k=1:rows
        u = yop.ast_control([name , '_', num2str(k)], weight(k), offset(k), deg(k));
        controls{k,1} = u;
        du = u.der;
        for n=1:deg
            controls{k,n+1} = du;
            du = du.der;
        end
    end
    varargout = cell(1, 1+deg);
    for k=1:deg+1
        varargout{k} = vertcat(controls{:,k});
    end
end

end