function u = control(varargin)

ip = inputParser();
ip.FunctionName = "yop.control";
ip.addOptional('size', [1, 1]);
ip.addParameter('name', 'u');
ip.addParameter('scaling', 1);
ip.addParameter('offset', 0);
ip.addParameter('int', 0);
ip.parse(varargin{:});

sz = ip.Results.size;
nu = prod(sz);
name = ip.Results.name;
scaling = ones(nu,1) .* ip.Results.scaling(:);
offset = ones(nu,1) .* ip.Results.offset(:);
int = ones(nu,1) .* ip.Results.int(:);

if any(int < 0)
    error(yop.error.invalid_control_degree());
end

if isequal(sz, [1, 1])
    u = yop.ast_control(name, scaling, offset, int);
else
   controls = cell(nu, 1);
    for k=1:nu
        controls{k} = yop.ast_control([name , '_', num2str(k)], ...
            scaling(k), offset(k), int(k));
    end
    u = reshape(vertcat(controls{:}), sz);
end

end