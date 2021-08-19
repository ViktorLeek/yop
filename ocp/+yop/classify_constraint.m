function c = classify_constraint(u)
% u - unclassified
% c - classified


% if isa(u.lhs, 'yop.ast_timepoint')
%     % Can access 'timepoint' property
%     if isa(u.timepoint, 'yop.ast_independent_initial')
%         return;
%     elseif isa(u.timepoint, 'yop.ast_independent_final')
%         return;
%     else
%         % inequality
%     end
% else
%
% end

switch class(u)
    case {'yop.ast_lt', 'yop.ast_le'}
        
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v <= 2
            c = yop.box_upper(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 <= v
            c = yop.box_lower(u.rhs, u.lhs);
            
        elseif isa(u.lhs, 'yop.ast_timepoint') && isnumeric(u.rhs)
            
            if isa(u.lhs.timpoint, 'yop.ast_independent_initial') && ...
               isa_variable(u.lhs.expr)
                % v(t0) < 2
                c  = yop.box_initial_upper(u.lhs.expr, u.rhs);
                
            else
                % Inequality constraint
            end
            
        else
            % Inequality constraint
            % lhs < rhs <==> lhs-rhs < 0
            c = yop.inequality_constraint(u.lhs-u.rhs);
            
        end
        
    case {'yop.ast_gt', 'yop.ast_ge'}
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v > 2
            c = yop.box_lower(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 > v
            c = yop.box_upper(u.rhs, u.lhs);
            
        else
            % Inequality constraint
            % lhs > rhs <==> rhs < lhs <==> rhs-lhs < 0
            c = yop.inequality_constraint(u.rhs-u.lhs);
            
        end
        
    case 'yop.ast_eq'
        if isa_variable(u.lhs) && isnumeric(u.rhs)
            % v == 2
            c = yop.box_equality(u.lhs, u.rhs);
            
        elseif isnumeric(u.lhs) && isa_variable(u.rhs)
            % 2 <= v
        end
        
    case 'yop.ast_ne'
    otherwise
        error('[yop] Error: Constraint type not regocnized');
end

end