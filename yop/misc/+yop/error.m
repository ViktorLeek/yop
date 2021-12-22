classdef error < handle
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
        
        function msg = failed_to_parse_ivp_eq()
            msg = [yop.error.msg_start(), ...
                'Failed to parse simulation equation', ...
                yop.error.msg_stop()];
        end
        
        function msg = initial_value_missing(var)
            msg = [yop.error.msg_start(), ...
                'No initial value for variable ' var.ast.name, ...
                'was given', yop.error.msg_stop()];
        end
        
        function msg = final_value_missing(var)
            msg = [yop.error.msg_start(), ...
                'No final value for variable ' var.ast.name, ...
                'was given', yop.error.msg_stop()];
        end
        
        function msg = start_of_integration_missing()
            msg = [yop.error.msg_start(), ...
                'Start of integration (t0) not specified', ...
                yop.error.msg_stop()];
        end
        
        function msg = end_of_integration_missing()
            msg = [yop.error.msg_start(), ...
                'End of integration (tf) not specified', ...
                yop.error.msg_stop()];
        end
        
        function msg = ivp_solver_not_recognized()
            msg = [yop.error.msg_start(), ...
                'IVP (simulation) solver not recognized', ...
                yop.error.msg_stop()];
        end
        
        function msg = ode_solver_for_dae_problem(solver)
            msg = [yop.error.msg_start(), ...
                'The problem cannot be solved using ', solver, ...
                ' because it is a DAE and ', solver, ' only hanldes ODEs', ...
                yop.error.msg_stop()];
        end
        
        function msg = ivp_grid_missing()
            msg = [yop.error.msg_start(), ...
                'No grid or sample points given for the simulation ', ...
                'problem', ...
                yop.error.msg_stop()];
        end
        
        function msg = yopvar_not_recognized(name)
            msg = [yop.error.msg_start(), ...
                'The ''yopvar'' ''', name, ''' could not be parsed.', ...
                newline, '''t'' is the independent variable, ''t0'' is its ', ...
                'initial value and ''tf'' is its finale value. States ', ...'
                'begin with ''x'', algebraics with ''z'', controls with ', ...
                '''u'', and parameters with ''p''' ...
                yop.error.msg_stop()];
        end
        
        
    end
end