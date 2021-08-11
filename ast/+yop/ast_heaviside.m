classdef ast_heaviside < yop.ast_functioncall
    properties (Constant)
        name = 'heaviside'
    end
    methods
        function obj = ast_heaviside(varargin)
            obj@yop.ast_functioncall(varargin{:});
        end
    end
end