function pc = sort_nonbox(con)

pc = yop.nbc_data();

for k=1:length(con)
    pck = con{k};
    switch class(pck)
        
        case 'yop.ast_eq'
            der_left = isa_der(pck.lhs);
            if all(der_left) % should add: all(isa_variable(pck.lhs))
                % der(v) == expr (assumed not der(v1) == der(v2))
                pc.add_ode(pck.lhs, pck.rhs);
                continue;
            elseif all(der_left == false) 
                % expr == ??, need to determine rhs
                r = pck; % remaining relations
            else
                % (der(x), expr) == (expr, ??)
                % mixed der, expr on lhs, rhs need to be determined for
                % expr part.
                tmp = yop.get_subrelation(pck, der_left);
                pc.add_ode(tmp.lhs, tmp.rhs);
                r = yop.ast_eq(yop.get_subrelation(pck, ~der_left));
            end
            
            der_right = isa_der(r.rhs);
            if all(der_right)
                pc.add_ode(r.rhs, r.lhs);
            elseif all(der_right == false)
                pc.add_eq(yop.ast_eq(r.lhs-r.rhs, 0));
            else
                tmp = yop.get_subrelation(r, der_right);
                pc.add_ode(tmp.rhs, tmp.lhs);
                pc.add_eq(yop.get_subrelation(r, ~der_right));
            end
            
        case {'yop.ast_le', 'yop.ast_lt'}
            c = yop.ast_le(pck.lhs - pck.rhs, 0);
            pc.add_ieq(c);
            
        case {'yop.ast_ge', 'yop.ast_gt'}
            c = yop.ast_le(pck.rhs - pck.lhs, 0);
            pc.add_ieq(c);
    end
end

end