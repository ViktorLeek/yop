function varargout = control(varargin)

ip = inputParser();
ip.FunctionName = "yop.control";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'u');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.addParameter('deg', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
nu = prod(sz);
name = ip.Results.name;
scaling = ones(nu,1) .* ip.Results.scaling(:);
offset = ones(nu,1) .* ip.Results.offset(:);
deg = ip.Results.deg;

if any(deg < 0)
    error(yop.error.invalid_control_degree());
end

if isequal(sz, [1, 1])
    u = yop.ast_control(name, scaling, offset, deg);
    varargout = {u};
    du = u.m_du;
    while ~isempty(du)
        varargout{end+1} = du;
        du = du.m_du;
    end
else
   controls = {}; 
    for k=1:nu
        u = yop.ast_control([name , '_', num2str(k)], scaling(k), offset(k), deg);
        controls{k,1} = u;
        du = u.m_du;
        for n=1:deg
            controls{k,n+1} = du;
            du = du.m_du;
        end
    end
    varargout = cell(1, 1+deg);
    for k=1:deg+1
        varargout{k} = reshape(vertcat(controls{:,k}), sz);
    end
end

end