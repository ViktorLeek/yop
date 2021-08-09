classdef subsref < yop.node
    
    properties
        arg
        subs
    end
    
    methods
        function obj = subsref(arg, subs)
            obj.arg = arg;
            obj.subs = subs;
            obj.dim = size( subsref( ones(size(arg)), subs ) );
        end
    end
    
    methods % Printing
        
        function print(obj)
            a = [];
            for k=1:length(obj.subs.subs)
                a = [a, 'a', num2str(k), ', '];
            end
            
            fprintf(['subsref(', a(1:end-2), ')\n']);
            
            % Here it is assumed that subs is vector of numerics, therefore
            % only the arg part is printed.
            last_child(obj);
            yop.print(obj.arg);
            end_child(obj);
        end
        
    end
end