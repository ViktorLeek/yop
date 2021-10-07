function [tps,ints,ders] = param_special_nodes(special_nodes, n_tp, n_int, N, tau, dt, t0, tf, t, x, u, p)
tps = [];
ints = [];
ders = [];
for node = special_nodes
    tmp_tp  = [tps;  zeros(n_tp  - length(tps), 1)];
    tmp_int = [ints; zeros(n_int - length(ints), 1)];
    switch node.type
        case yop.ocp_expr.tp
            tp = yop.parameterize_timepoint(node,t0,tf,t,x,u,p,tmp_tp,tmp_int,ders);
            tps = [tps; tp(:)];
            
        case yop.ocp_expr.int
            int = yop.parameterize_integral(node,N,tau,dt,t,x,u,p,tmp_tp,tmp_int,ders);
            ints = [ints; int(:)];
            
        case yop.ocp_expr.der
            error(yop.msg.not_implemented);
            warning(['When implementing dont forget to ' ...
                'multiply derivative with step length']);
        otherwise
            error(yop.msg.unexpected_error);
    end
end
end