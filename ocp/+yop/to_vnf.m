function vnf = to_vnf(svhsrf, vnf)
% To variable-numeric form
% Last step before it can be distinguhised which constraints are box and
% which are other constraints

if nargin == 1
    vnf = yop.vnf_data();
end

for k=1:length(svhsrf.ve)
    rk = svhsrf.ve{k};
    cnstrctr = yop.get_constructor(rk);
    
    nums = isa_numeric(rk.rhs);
    
    if all(nums)
        vnf.add_vn(rk);
        
    elseif all(~nums)
        vnf.add_ve(rk);
        
    else
        vnf.add_vn(cnstrctr(rk.lhs( nums), rk.rhs( nums)));
        vnf.add_ve(cnstrctr(rk.lhs(~nums), rk.rhs(~nums)));
    end 
end

for k=1:length(svhsrf.ev)
    rk = svhsrf.ev{k};
    cnstrctr = yop.get_constructor(rk);
    
    nums = isa_numeric(rk.lhs);
    
    if all(nums)
        vnf.add_nv(rk);
        
    elseif all(~nums)
        vnf.add_ev(rk);
        
    else
        vnf.add_nv(cnstrctr(rk.lhs( nums), rk.rhs( nums)));
        vnf.add_ev(cnstrctr(rk.lhs(~nums), rk.rhs(~nums)));
    end 
end

end