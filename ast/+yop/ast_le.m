classdef ast_le < yop.ast_relation
    
    properties (Constant)
        name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
            %obj.dim = le(ones(size(lhs)), ones(size(rhs)));
        end
        
        function value = evaluate(obj)
            value = le(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = le(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
        
        function c = parse(obj)
            
            % obj is assumed to be on single relation form. To ensure that,
            % call the function 'yop.to_srf' which gives all relevenat
            % subrelations on srf form. obj is also assumed to be on
            
            
            % Identify box constraint:  v <= num
            if all(isa_variable(obj.lhs))
                if isnumeric(obj.rhs)
                    c = yop.box_upper(obj.lhs, obj.rhs);
                    return;
                end
            end
            
            % Identify box constraint:  num <= v
            if isnumeric(obj.lhs)
                if all(isa_variable(obj.rhs))
                    c = yop.box_lower(obj.rhs, obj.lhs);
                    return;
                end
            end
            
            % Identify box constratint: var(tp) <= num, tp in {t, t0, tf}
            if isa(obj.lhs, 'yop.ast_timepoint')       % some(tp) <= some
                if isnumeric(obj.rhs)                  % some(tp) <= num
                    if all(isa_variable(obj.lhs.expr)) % var(tp)  <= num
                        
                        if isa(obj.lhs.timepoint, 'yop.ast_independent')
                            c = yop.box_upper(obj.lhs.expr, obj.rhs);
                            return;
                        end
                        
                        if isa(obj.lhs.timepoint, 'yop.ast_independent_initial')
                            c  = yop.box_initial_upper(obj.lhs.expr, obj.rhs);
                            return;
                        end
                        
                        if isa(obj.lhs.timepoint, 'yop.ast_independent_final')
                            c = yop.box_final_upper(obj.lhs.expr, obj.rhs);
                            return;
                        end
                    end
                end
            end
            
            % Identify box constratint: num <= var(tp), tp in {t, t0, tf}
            if isa(obj.rhs, 'yop.ast_timepoint')
                if isnumeric(obj.lhs)
                    if all(isa_variable(obj.rhs.expr))
                        
                        if isa(obj.rhs.timepoint, 'yop.ast_independent')
                            c = yop.box_lower(obj.rhs.expr, obj.lhs);
                            return;
                        end
                        
                        if isa(obj.rhs.timepoint, 'yop.ast_independent_initial')
                            c = yop.box_initial_lower(obj.rhs.expr, obj.lhs);
                            return;
                        end
                        
                        if isa(obj.rhs.timepoint, 'yop.ast_independent_final')
                            c = yop.box_final_lower(obj.rhs.expr, obj.lhs);
                            return;
                        end
                    end
                end
            end
            
            % Arbitrary inequality constraint: expr <= 0
            c = yop.inequality_constraint(obj.lhs-obj.rhs);
            
        end
    end
end