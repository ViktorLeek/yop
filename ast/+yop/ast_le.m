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
        
        function rels = split_vars(srf)
            % Transformations before: srf, split vars and exprs
            
            if all(isa_variable(srf.lhs)) && isnumeric(srf.rhs)
                
                % In order to split the variables it is first necessary to
                % deterimine the variables that reaches the definition of
                % srf.lhs
                rv = reaching_variables(srf.lhs);
                
                % Based on the reaching variables the expression is divided
                % into subexpressions.
                rels = {};
                for k=1:length(rv)
                    if ~isempty(rv(k).reaching)
                        idx = rv(k).index;
                        if isscalar(srf.rhs)
                            node = yop.ast_le(srf.lhs(idx), srf.rhs);
                            
                        else
                            node = yop.ast_le(srf.lhs(idx), srf.rhs(idx));
                            
                        end
                        rels = {rels{:}, node};  
                    end
                end
                assert(~isempty(rels), ...
                    ['[Yop] Unexpected error: Some variables are', ...
                    'expected to reach the expression. Check if the', ...
                    'proper transforms has been applied first.'])
                
            elseif isnumeric(srf.lhs) && all(isa_variable(srf.rhs))
                
                rv = reaching_variables(srf.rhs);
                
                % Based on the reaching variables the expression is divided
                % into subexpressions.
                rels = {};
                for k=1:length(rv)
                    if ~isempty(rv(k).reaching)
                        idx = rv(k).index;
                        if isscalar(srf.lhs)
                            node = yop.ast_le(srf.lhs, srf.rhs(idx));
                            
                        else
                            node = yop.ast_le(srf.lhs(idx), srf.rhs(idx));
                            
                        end
                        rels = {rels{:}, node};  
                    end
                end
                assert(~isempty(rels), ...
                    ['[Yop] Unexpected error: Some variables are', ...
                    'expected to reach the expression. Check if the', ...
                    'proper transforms has been applied first.'])
                
            else
                % If it is neither var <= num nor num <= var then this
                % transformation is unnecessary because the relation is not
                % a box constraint.
                rels = {srf};
            end
        end
        
        function rels = split_vars_and_exprs(srf)
            % rels = split_vars_and_exprs(srf)
            % split the relation 'srf' into relations that do not mix
            % expressions and variables on the same side. It is either:
            %   1) var  <= var
            %   2) var  <= expr
            %   3) expr <= var
            %   4) expr <= expr
            %
            % Never: [expr; var] <= some_value
            
            % First process lhs
            vars = isa_variable(srf.lhs);
            
            if all(vars) || all(~vars)
                % No need to split if all are or none are variables
                rels_tmp = {srf};
                
            elseif isscalar(srf.rhs)
                % Cannot access subindices of scalars
                rels_tmp = {...
                    yop.ast_le(srf.lhs( vars), srf.rhs), ...
                    yop.ast_le(srf.lhs(~vars), srf.rhs)
                    };
            else
                rels_tmp = {...
                    yop.ast_le(srf.lhs( vars), srf.rhs( vars)), ...
                    yop.ast_le(srf.lhs(~vars), srf.rhs(~vars))
                    };
            end
            % At this point, there should be no mix of variables and
            % expressions on lhs.
            
            % Process rhs
            rels = {};
            for k=1:length(rels_tmp)
                rk = rels_tmp{k};
                vars = isa_variable(rk.rhs);
                
                if all(vars) || all(~vars)
                    % No need to split if all are or none are variables
                    rels = {rels{:}, rk};
                    
                elseif isscalar(rk.lhs)
                    % Cannot access subindices of scalars
                    rels = {rels{:}, ...
                        yop.ast_le(srf.lhs, srf.rhs( vars)), ...
                        yop.ast_le(srf.lhs, srf.rhs(~vars))
                        };
                else
                    rels = {rels{:}, ...
                        yop.ast_le(srf.lhs( vars), srf.rhs( vars)), ...
                        yop.ast_le(srf.lhs(~vars), srf.rhs(~vars))
                        };
                end
            end
            
            for k=1:length(rels)
                rk = rels{k};
                vars = isa_variable(rk.lhs);
                assert(all(vars) || all(~vars));
                vars = isa_variable(rk.rhs);
                assert(all(vars) || all(~vars));
            end
            
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