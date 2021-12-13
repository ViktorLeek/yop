classdef error < handle
    properties
        msg
    end
    methods
        function obj = error(msg)
            obj.msg = msg;
        end
    end
    
    methods (Static)
        function m = msg_start()
            m = '[Yop] ';
        end
        
        function m = msg_stop()
            m = '.';
        end
        
        function msg = timevarying_objective(err_vars)
            names = '';
            for k=1:length(err_vars)
                names = [names, '''', err_vars{k}.name, '''', ', '];
            end
            names = names(1:end-2);
            
            msg = [yop.error.msg_start(), ...
                'Objective function error. ', newline, ...
                'Variable(s): ', names, ' is not part of an integral or ' ...
                'timepoint and can therefore not be evaluted in an ' ...
                'objective function context. Consider evaluating ' ...
                'inside an integral or at a timepoint instead. ' ...
                'Valid objective functions are '...
                'functions of timepoints, integrals and parameters', ...
                yop.error.msg_stop()];
        end
        
        function msg = multiple_independent_variables()
            msg = [yop.error.msg_start(), ...
                'Detected multiple independent variables', ...
                yop.error.msg_stop()];
        end
        
        function multiple_independent_initial()
            msg = [yop.error.msg_start(), ...
                'Detected multiple initial values for the ', ...
                'independent variable', ...
                yop.error.msg_stop()];
            yop.errors.report(yop.error(msg));
        end
        
        function msg = multiple_independent_final()
            msg = [yop.error.msg_start(), ...
                'Detected multiple final values for the ', ...
                'independent variable', ...
                yop.error.msg_stop()];
        end
        
        function msg = incompatible_constraint_size()
            msg = [yop.error.msg_start(), ...
                'Incompatible size for constraint lhs and rhs', ...
                yop.error.msg_stop()];
        end
        
        function msg = unknown_constraint()
            msg = [yop.error.msg_start(), ...
                'Unknown constrain type', ...
                yop.error.msg_stop()];
        end
        
        function msg = failed_to_find_variable(id)
            msg = [yop.error.msg_start(), ...
                'Failed to find variable ', ...
                'with id [', num2str(id), ']', ...
                yop.error.msg_stop()];
        end
        
        function msg = missing_state_derivative(err_vars)
            names = '';
            for k=1:length(err_vars)
                names = [names, '''', err_vars{k}.name, '''', ', '];
            end
            names = names(1:end-2);
            
            msg = [yop.error.msg_start(), ...
                'Differential equation error. ', newline, ...
                'State(s): ', names, ' is(are) not given a ', ...
                'differential equation' ...
                yop.error.msg_stop()];
        end
        
        function msg = failed_to_parse_timed_expression()
            msg = [yop.error.msg_start(), ...
                'Failed to parse timepoint or interval', ...
                yop.error.msg_stop()];
        end
        
    end
end