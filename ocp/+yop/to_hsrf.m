function hsrf = to_hsrf(srf_cell, hsrf)
% Transform representation to homogeneous single relation form
% split the srf into relations that do not mix
% expressions and variables on the same side. It is either:
%   1) vv: var  <= var
%   2) ve: var  <= expr
%   3) ev: expr <= var
%   4) ee: expr <= expr
%   :) u is for unknown
% Never: [expr; var] <= some_value
%
% srf - A single relation on srf form
% rel - function handle to the

if nargin == 1
    hsrf = yop.hsrf_data();
else
    hsrf.clear_unknown();
end

for n=1:length(srf_cell)
    rn = srf_cell{n};
    
    % Function handle to the class contructor
    cnstrctr = yop.get_constructor(rn);
    
    % First process lhs
    vars = isa_variable(rn.lhs);
    
    if all(vars)
        % No need to split if all are variables
        hsrf.add_vu(rn);
        
    elseif all(~vars)
        % No need to split if none are variables
        hsrf.add_eu(rn);
        
    elseif isscalar(rn.rhs)
        % Cannot access subindices of scalars
        hsrf.add_vu(cnstrctr(rn.lhs( vars), rn.rhs));
        hsrf.add_eu(cnstrctr(rn.lhs(~vars), rn.rhs));
    else
        hsrf.add_vu(cnstrctr(rn.lhs( vars), rn.rhs( vars)));
        hsrf.add_eu(cnstrctr(rn.lhs(~vars), rn.rhs(~vars)));
    end
    
    
    for k=1:length(hsrf.vu)
        rk = hsrf.vu{k};
        vars = isa_variable(rk.rhs);
        
        if all(vars)
            hsrf.add_vv(rk);
            
        elseif all(~vars)
            hsrf.add_ve(rk);
            
        elseif isscalar(rk.lhs)
            hsrf.add_vv(cnstrctr(rk.lhs, rk.rhs( vars)));
            hsrf.add_ve(cnstrctr(rk.lhs, rk.rhs(~vars)));
        else
            hsrf.add_vv(cnstrctr(rk.lhs( vars), rk.rhs( vars)));
            hsrf.add_ve(cnstrctr(rk.lhs(~vars), rk.rhs(~vars)));
        end
    end
    
    for k=1:length(hsrf.eu)
        rk = hsrf.eu{k};
        vars = isa_variable(rk.rhs);
        
        if all(vars)
            hsrf.add_ev(rk);
            
        elseif all(~vars)
            hsrf.add_ee(rk);
            
        elseif isscalar(rk.lhs)
            hsrf.add_ev(cnstrctr(rk.lhs, rk.rhs( vars)));
            hsrf.add_ee(cnstrctr(rk.lhs, rk.rhs(~vars)));
        else
            hsrf.add_ev(cnstrctr(rk.lhs( vars), rk.rhs( vars)));
            hsrf.add_ee(cnstrctr(rk.lhs(~vars), rk.rhs(~vars)));
        end
    end
    
    hsrf.clear_unknown();
end
end