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
        dtp.add_vn_t( yop.get_subrelation(rk, ~is_tp) );
    end
    
    if ~all(~t0) % Do not add if t0 is all false
        dtp.add_vn_t0( yop.get_subrelation(rk, t0) );
    end
    
    if ~all(~tf) 
        dtp.add_vn_tf( yop.get_subrelation(rk, tf) );
    end
    
    if ~all(~tx)
        % Numerical timepoint are not box constraints (yet, at least) so
        % they are moved to the general category of expression-variable
        % relation
        dtp.add_ve( yop.get_subrelation(rk, tx) );
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
        dtp.add_nv_t( yop.get_subrelation(rk, ~is_tp) );
    end
    
    if ~all(t0==false)
        dtp.add_nv_t0( yop.get_subrelation(rk, t0) );
    end
    
    if ~all(tf==false) 
        dtp.add_nv_tf( yop.get_subrelation(rk, tf) );
    end
    
    if ~all(tx==false)
        % Numerical timepoint are not box constraints (yet, at least) so
        % they are moved to the general category of expression-variable
        % relation
        dtp.add_ev( yop.get_subrelation(rk, tx) );
    end
    
    assert(all(~is_tp | t0 | tf | tx), '[Yop] Unexpected error.');
end

end