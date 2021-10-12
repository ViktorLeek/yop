function varargout = control(varargin)

ip = inputParser();
ip.FunctionName = "yop.control";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'u');
ip.addParameter('pw', 'constant');
ip.parse(varargin{:});

sz = ip.Results.size;
name = ip.Results.name;
pw = ip.Results.pw;

if isscalar(sz)
    sz = [sz, 1];
end

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

u = yop.ast_control(name, sz(1), sz(2), pw);

varargout = {u};
du = u.der;
while ~isempty(du)
    varargout{end+1} = du;
    du = du.der;
end

end