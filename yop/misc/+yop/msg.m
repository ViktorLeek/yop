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
        
        function m = ivp_relation()
            m = [yop.msg.start, ['IVPs can only have equality ' ,...
                'relations'], yop.msg.stop];
        end
        
        function m = ivp_ub_differ()
            m = [yop.msg.start, ['The initial value of the ', ...
                'independent variable (time) is not exactly determined'],...
                yop.msg.stop];
        end
        
        function m = ivp_no_start_time()
            m = [yop.msg.start, ['The initial value of the ', ...
                'independent variable (time) is not set'],...
                yop.msg.stop];
        end
        
        function m = ivp_t0_err()
            m = [yop.msg.start, ['The timepoint set for the variable ', ...
                'does not correspond to the initial timepoint'],...
                yop.msg.stop];
        end
        
    end
end