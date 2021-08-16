function error(msg, varargin)
err_msg = ['[yop] Error: ', msg];
error(err_msg, varargin{:});
end