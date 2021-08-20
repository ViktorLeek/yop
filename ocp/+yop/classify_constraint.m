function c = classify_constraint(u)
% classify_constraint(u) - Classifies the constraint u. u is assumed to be
% in single relation form (srf).
% u - unclassified
% c - classified

switch class(u)
    case {'yop.ast_lt', 'yop.ast_le'}
        % Will be undefined if der(x) < expr
        
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v < 2
            c = yop.box_upper(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 < v
            c = yop.box_lower(u.rhs, u.lhs);
            
        elseif isa(u.lhs, 'yop.ast_timepoint') && isnumeric(u.rhs)
            % expr(tp) < 2
            
            if isa(u.lhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.lhs.expr)
                % v(t) < 2
                c = yop.box_upper(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.lhs.expr)
                % v(t0) < 2
                c  = yop.box_initial_upper(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.lhs.expr)
                % v(tf) < 2
                c = yop.box_final_upper(u.lhs.expr, u.rhs); 

            else
                % Inequality constraint
                % expr(tp) < 2 <==> expr(tp) - 2 < 0
                c = yop.inequality_constraint(u.lhs-u.rhs);
                
            end
            
        elseif isnumeric(u.lhs) && isa(u.rhs, 'yop.ast_timepoint')
            % 2 < expr(tp)
            
            if isa(u.rhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.rhs.expr)
                % 2 < v(t)
                c = yop.box_lower(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.rhs.expr)
                % 2 < v(t0) <==> v(t0) > 2
                c = yop.box_initial_lower(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.rhs.expr)
                % 2 < v(tf) <==> v(tf) > 2
                c = yop.box_final_lower(u.rhs.expr, u.lhs);
                
            else 
                % Inequality constraint
                % 2 < expr(tp) <==> 2 - expr(tp) < 0
                c = yop.inequality_constraint(u.lhs - u.rhs);
            end
            
            
        else
            % Inequality constraint
            % lhs < rhs <==> lhs-rhs < 0
            c = yop.inequality_constraint(u.lhs-u.rhs);
            
        end
        
    case {'yop.ast_gt', 'yop.ast_ge'}
        % Will be undefined if der(x) > expr
        
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v > 2
            c = yop.box_lower(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 > v
            c = yop.box_upper(u.rhs, u.lhs);
            
        elseif isa(u.lhs, 'yop.ast_timepoint') && isnumeric(u.rhs)
            % expr(tp) > 2
            
            if isa(u.lhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.lhs.expr) 
                % v(t) > 2
                c = yop.box_lower(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.lhs.expr)
                % v(t0) > 2
                c  = yop.box_initial_lower(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.lhs.expr)
                % v(tf) > 2
                c = yop.box_final_lower(u.lhs.expr, u.rhs); 

            else
                % Inequality constraint
                % expr(tp) > 2 <==> 2 < expr(tp) <==> 2 - expr(tp) < 0
                c = yop.inequality_constraint(u.rhs-u.lhs);
            end
            
        elseif isnumeric(u.lhs) && isa(u.rhs, 'yop.ast_timepoint')
            % 2 > expr(tp)
            
            if isa(u.rhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.rhs.expr)
                % 2 > v(t) <==> v(t) < 2
                c = yop.box_upper(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.rhs.expr)
                % 2 > v(t0) <==> v(t0) < 2
                c = yop.box_initial_upper(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.rhs.expr)
                % 2 > v(tf) <==> v(tf) < 2
                c = yop.box_final_upper(u.rhs.expr, u.lhs);
                
            else 
                % Inequality constraint
                % 2 > expr(tp) <==> expr(tp) < 2 <==> expr(tp) - 2 < 0
                c = yop.inequality_constraint(u.rhs - u.lhs);
            end
            
        else
            % Inequality constraint
            % lhs > rhs <==> rhs < lhs <==> rhs-lhs < 0
            c = yop.inequality_constraint(u.rhs-u.lhs);
            
        end
        
    case 'yop.ast_eq'
        % der, alg
        
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v == 2
            c = yop.box_equality(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 <= v
            c = yop.box_equality(u.rhs, u.lhs);
            
        elseif isa(u.lhs, 'yop.ast_timepoint') && isnumeric(u.rhs)
            % expr(tp) = 2
            
            if isa(u.lhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.lhs.expr) 
                % v(t) = 2
                c = yop.box_equality(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.lhs.expr)
                % v(t0) = 2
                c  = yop.box_initial_equality(u.lhs.expr, u.rhs);
                
            elseif isa(u.lhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.lhs.expr)
                % v(tf) = 2
                c = yop.box_final_equality(u.lhs.expr, u.rhs); 

            else
                % Inequality constraint
                % expr(tp) = 2 <==> expr(tp) - 2 = 0
                c = yop.equality_constraint(u.lhs-u.rhs);
            end
            
        elseif isnumeric(u.lhs) && isa(u.rhs, 'yop.ast_timepoint')
            % 2 = expr(tp)
            
            if isa(u.rhs.timepoint, 'yop.ast_independent') && ...
                    isa_variable(u.rhs.expr)
                % 2 = v(t)
                c = yop.box_equality(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_initial') && ...
                    isa_variable(u.rhs.expr)
                % 2 = v(t0)
                c = yop.box_initial_equality(u.rhs.expr, u.lhs);
                
            elseif isa(u.rhs.timepoint, 'yop.ast_independent_final') && ...
                    isa_variable(u.rhs.expr)
                % 2 = v(tf)
                c = yop.box_final_equality(u.rhs.expr, u.lhs);
                
            else 
                % Inequality constraint
                % 2 = expr(tp) <==> 2 - expr(tp) = 0
                c = yop.equality_constraint(u.lhs - u.rhs);
            end
            
        elseif isa(u.lhs, 'yop.ast_der') 
            c = yop.differential_constraint(u.lhs.var, u.rhs);
            
        elseif isa(u.rhs, 'yop.ast_der')
            c = yop.differential_constraint(u.rhs.var, u.lhs);
            
        elseif isa(u.lhs, 'yop.ast_alg')
            c = yop.algebraic_contraint(u.lhs);
            
        elseif isa(u.rhs, 'yop.ast_alg')
            c = yop.algebraic_contraint(u.rhs);
            
        else
            c = yop.equality_constraint(u.lhs - u.rhs);
            
        end
        
    case 'yop.ast_ne'
        error('[yop] Error: Yop does not allow not equal as a constraint');
        
    otherwise
        error('[yop] Error: Constraint type not recognized');
end

end