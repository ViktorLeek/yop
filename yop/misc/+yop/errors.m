classdef errors < handle
    properties
        errs = yop.error.empty(1,0)
    end
    methods
        function obj = errors()
            persistent OBJ
            if isempty(OBJ)
                OBJ = obj;
            else
                obj = OBJ;
            end
        end
    end
    
    methods (Static)
        function clear()
            e = yop.errors();
            e.errs = yop.error.empty(1,0);
        end
        
        function report(error)
            e = yop.errors();
            e.errs(end+1) = error;  
        end
        
        function bool = empty()
            e = yop.errors();
            bool = isempty(e.errs);
        end
        
        function msg = get()
            % Get and forget (implicit clear)
            % The reason for the implicit clear is that it is assumed that
            % this call is made just before issuing an error, after which
            % it is desired that all old errors are forgotten.
            msg = [];
            for e = yop.errors().errs
                msg = [msg, e.msg, newline];
            end
            yop.errors.clear();
        end
    end
end