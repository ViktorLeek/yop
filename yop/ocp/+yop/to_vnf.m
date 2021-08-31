function vnf = to_vnf(hsrf, vnf)
% To variable-numeric form
% Last step before it can be distinguhised which constraints are box and
% which are other constraints

if nargin == 1
    vnf = yop.vnf_data();
end

% Even though not all relations are affected by the transformation, the
% relations still represent the complete constraint set, so all relations
% that are not affected are copied.
vnf.add_vv(hsrf.vv);
vnf.add_ee(hsrf.ee);

for k=1:length(hsrf.ve)
    rk = hsrf.ve{k};
    
    nums = isa_numeric(rk.rhs);
    
    if all(nums)
        vnf.add_vn(rk);
        
    elseif all(~nums)
        vnf.add_ve(rk);
        
    else
        vnf.add_vn( yop.get_subrelation(rk,  nums) );
        vnf.add_ve( yop.get_subrelation(rk, ~nums) );
    end 
    
end

for k=1:length(hsrf.ev)
    rk = hsrf.ev{k};
    
    nums = isa_numeric(rk.lhs);
    
    if all(nums)
        vnf.add_nv(rk);
        
    elseif all(~nums)
        vnf.add_ev(rk);
        
    else
        vnf.add_nv( yop.get_subrelation(rk,  nums) );
        vnf.add_ev( yop.get_subrelation(rk, ~nums) );
    end 
end

end

% function sr = get_subrelation(relation, idx)
% % Since indices might be scaled and variables can be scalars it is
% % necessary to test if it is possible to take the subindices of the
% % relations.
% 
% if isscalar(relation.rhs)
%     rhs = relation.rhs;
% else
%     rhs = relation.rhs(idx);
% end
% 
% if isscalar(relation.lhs)
%     lhs = relation.lhs;
% else
%     lhs = relation.lhs(idx);
% end
% 
% f  = get_constructor(relation);
% sr = f(lhs, rhs);
% 
% end