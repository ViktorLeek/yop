function dtp = to_dtp(vnf)
% Transform the vnf form to relations with distinct timepoint.
% Only concerns the sets vn and nv

dtp = yop.dtp_data();

% preserve all old analysis information
% Should be dtp_data(vnf);
dtp.add_vv(vnf.vv);
dtp.add_ve(vnf.ve);
dtp.add_ev(vnf.ev);
dtp.add_ee(vnf.ee);

for k=1:length(vnf.vn)
    rk = vnf.vn{k}; % rk - relation k
    
    [is_tp, tp] = isa_timepoint(rk.lhs);
    
    % The elements to select
    t0 = tp==yop.initial_timepoint;
    tf = tp==yop.final_timepoint;
    tx = ~t0 & ~tf & is_tp; % Numerical timepoints
    
    if ~all(is_tp)
        dtp.add_vn_t( get_subrelation(rk, ~is_tp) );
    end
    
    if ~all(~t0) % Do not add if t0 is all false
        dtp.add_vn_t0( get_subrelation(rk, t0) );
    end
    
    if ~all(~tf) 
        dtp.add_vn_tf( get_subrelation(rk, tf) );
    end
    
    if ~all(~tx)
        % Numerical timepoint are not box constraints (yet, at least) so
        % they are moved to the general category of expression-variable
        % relation
        dtp.add_ve( get_subrelation(rk, tx) );
    end
    
    assert(all(~is_tp | t0 | tf | tx), '[Yop] Unexpected error.');
end

for k=1:length(vnf.nv)
    rk = vnf.nv{k}; % rk - relation k
    
    [is_tp, tp] = isa_timepoint(rk.rhs);
    
    % The elements to select
    t0 = tp==yop.initial_timepoint;
    tf = tp==yop.final_timepoint;
    tx = ~t0 & ~tf & is_tp; % Numerical timepoints
    
    if ~all(is_tp)
        dtp.add_nv_t( get_subrelation(rk, ~is_tp) );
    end
    
    if ~all(~t0) % Do not add if t0 is all false
        dtp.add_nv_t0( get_subrelation(rk, t0) );
    end
    
    if ~all(~tf) 
        dtp.add_nv_tf( get_subrelation(rk, tf) );
    end
    
    if ~all(~tx)
        % Numerical timepoint are not box constraints (yet, at least) so
        % they are moved to the general category of expression-variable
        % relation
        dtp.add_ev( get_subrelation(rk, tx) );
    end
    
    assert(all(~is_tp | t0 | tf | tx), '[Yop] Unexpected error.');
end

end

function sr = get_subrelation(relation, idx)
% Since indices might be scaled and variables can be scalars it is
% necessary to test if it is possible to take the subindices of the
% relations.

if isscalar(relation.rhs)
    rhs = relation.rhs;
else
    rhs = relation.rhs(idx);
end

if isscalar(relation.lhs)
    lhs = relation.lhs;
else
    lhs = relation.lhs(idx);
end

f  = get_constructor(relation);
sr = f(lhs, rhs);

end