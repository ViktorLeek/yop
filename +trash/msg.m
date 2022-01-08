classdef msg < handle
   
    % m = [yop.msg.start, '', yop.msg.stop];
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
        
        function m = control_degree_error()
            m = [yop.msg.start, 'Invalid option for the control ' ...
                'parametrization. Set ''pw'' to either of ''constant'','...
                ' ''linear'', ''quadratic'', ''qubic'', or to a ' ...
                'non-negative number, 0 corresponds to ''constant'', ' ...
                '1 to ''linear'', and so on, with no upper limit', ...
                yop.msg.stop];
        end
        
        function m = ival_relation_error(clss)
            m = [yop.msg.start, 'The only relations allowd in a time ' ...
                'intervals are relations. You have ''',  clss, '''', ...
                yop.msg.stop];
        end
        
        function m = timed_expr_ill_expr()
            m = [yop.msg.start, 'Illegal expression for timepoint or ' ...
                'time interval detected', yop.msg.stop];
        end
        
        function m = timed_expr_relation_overflow()
            m = [yop.msg.start, 'Too many relations for a timepoint or '...
                'time interval', yop.msg.stop];
        end
        
        function m = illegal_timepoint()
            m = [yop.msg.start, 'Cannot parse the supplied timepoint', ...
                yop.msg.stop];
        end
        
        function m = ambig_ival()
            m = [yop.msg.start, 'Time interval is ambiguous', ...
                yop.msg.stop];
        end
        
        function m = err_objective(vars)
            names = '';
            for k=1:length(vars)
                names = [names, '''', vars{k}.name, '''', ', '];
            end
            names = names(1:end-2);
            
            m = [yop.msg.start, 'Objective function error. ', newline, ...
                'Variable(s): ', names, ' is not part of an integral or ' ...
                'timepoint and can therefore not be evaluted in an ' ...
                'objective function context. Consider evaluating ' ...
                'inside an integral or at a timepoint instead. ' ...
                'Valid objective functions can only '...
                'be functions of timepoints, integrals and parameters', ...
                yop.msg.stop];
        end
        
    end
end