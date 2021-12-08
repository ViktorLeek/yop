function varargout = control(varargin)

ip = inputParser();
ip.FunctionName = "yop.control";
ip.addOptional('rows', 1);
ip.addParameter('name', 'u');
ip.addParameter('pw', 'constant');
ip.parse(varargin{:});

rows = ip.Results.rows;
name = ip.Results.name;
pw = ip.Results.pw;

switch pw
    case 'constant'
        pw = 0;
    case 'linear'
        pw = 1;
    case 'quadratic'
        pw = 2;
    case 'cubic'
        pw = 3;
    otherwise
        if ~isnumeric(pw) || ~isscalar(pw) || any(pw < 0)
            error(yop.msg.control_degree_error);
        end
end

if rows == 1
    u = yop.ast_control(name, pw);
    varargout = {u};
    du = u.der;
    while ~isempty(du)
        varargout{end+1} = du;
        du = du.der;
    end
else
   controls = cell(rows, 1+pw); 
    for k=1:rows
        u = yop.ast_control([name , '_', num2str(k)], pw);
        controls{k,1} = u;
        du = u.der;
        for n=1:pw
            controls{k,n+1} = du;
            du = du.der;
        end
    end
    varargout = cell(1, 1+pw);
    for k=1:pw+1
        varargout{k} = vertcat(controls{:,k});
    end
end

end