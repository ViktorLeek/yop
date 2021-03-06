function varargout = vars(varargin)
ip = inputParser();
ip.FunctionName = "yop.vars";
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
    varargout{end+1} = yop.state(yop.IF(isempty(nx), [], [nx,1]));
end

if ~isempty(nz)
    varargout{end+1} = yop.algebraic(yop.IF(isempty(nz), [], [nz,1]));
end

if ~isempty(nu)
    varargout{end+1} = yop.control(yop.IF(isempty(nu), [], [nu,1]));
end

if ~isempty(np)
    varargout{end+1} = yop.parameter(yop.IF(isempty(np), [], [np,1]));
end

end