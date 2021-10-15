function varargout = ocp_variables(varargin)
ip = inputParser();
ip.FunctionName = "yop.ocp_variables";
ip.addParameter('nx', []);
ip.addParameter('nz', []);
ip.addParameter('nu', []);
ip.addParameter('np', []);
ip.parse(varargin{:});

nx  = ip.Results.nx;
nz  = ip.Results.nz;
nu  = ip.Results.nu;
np  = ip.Results.np;

t0 = yop.time0();
tf = yop.timef();
t = yop.time();
varargout = {t0, tf, t};

if ~isempty(nx)
    varargout{end+1} = yop.state(nx);
end

if ~isempty(nz)
    varargout{end+1} = yop.algebraic(nz);
end

if ~isempty(nu)
    varargout{end+1} = yop.control(nu);
end

if ~isempty(np)
    varargout{end+1} = yop.parameter(np);
end

end