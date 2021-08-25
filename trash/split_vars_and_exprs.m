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
            
            if all(vars) 
                % No need to split if all are variables
                rels_tmp = {srf};
                
            elseif all(~vars)
                % No need to split if none are variables
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
                    % No need to split if all or none are variables
                    rels = {rels{:}, rk};
                    
                elseif isscalar(rk.lhs)
                    % Cannot access subindices of scalars
                    rels = {rels{:}, ...
                        yop.ast_le(rk.lhs, rk.rhs( vars)), ...
                        yop.ast_le(rk.lhs, rk.rhs(~vars))
                        };
                else
                    rels = {rels{:}, ...
                        yop.ast_le(rk.lhs( vars), rk.rhs( vars)), ...
                        yop.ast_le(rk.lhs(~vars), rk.rhs(~vars))
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