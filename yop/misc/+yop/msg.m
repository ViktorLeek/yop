classdef msg < handle
   
    methods (Static)
        
        function m = start()
            m = '[Yop] ';
        end
        
        function m = stop()
            m = '.';
        end
        
        function m = default_der(name)
            m = [yop.msg.start, 'Some or all elemnts of state "', name, ...
                '" is not bound by a ODE. Setting default value: der(', ...
                name, ') == 0, for elements that are not ', ...
                'already set', yop.msg.stop];
        end
        
        function m = not_implemented()
            m=[yop.msg.start,'Using a feature that is not implemented ',...
                yop.msg.stop];
        end
        
        function m = unexpected_error()
            m = [yop.msg.start, 'Unexpected error', yop.msg.stop];
        end
        
        
    end
end