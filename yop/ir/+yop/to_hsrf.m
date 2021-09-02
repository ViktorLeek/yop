function hsrf = to_hsrf(relations, hsrf)
% Transform representation to homogeneous single relation form.
%   Split the srf into relations that do not mix expressions and variables 
%   on the same side. It is either:
%     1) vv: var  <= var
%     2) ve: var  <= expr
%     3) ev: expr <= var
%     4) ee: expr <= expr
%     :) u is for unknown
%   Never: [expr; var] <= some_value
% 
%   Furthermore, it also separates the variables so that ve and ev 
%   relations only contains variables with a unique ID

if nargin == 1
    hsrf = yop.hsrf_data();
else
    hsrf.clear_unknown();
end

for n=1:length(relations)
    %% Split relation into var-unknowns and expr-unknowns
    rn = relations{n};
    [var, ID] = isa_variable_helper(rn, 'lhs');
    
    % Add variable-unknown relation
    for u = yop.row_vec( unique(ID(var)) )
        hsrf.add_vu( yop.get_subrelation(rn, ID==u) );
    end
    
    % Add expression-unknown relation
    if ~all(var)
        hsrf.add_eu( yop.get_subrelation(rn, ~var) );
    end
    
    %% Split var-unknown into var-var and var-expr
    % Process all variable-unknown
    for k=1:length(hsrf.vu)
        ck = hsrf.vu{k};
        [var, ID] = isa_variable_helper(ck, 'rhs');
        for u = yop.row_vec(unique(ID(var)))
            hsrf.add_vv( yop.get_subrelation(ck, ID==u) );
        end
        
        if ~all(var)
            hsrf.add_ve( yop.get_subrelation(ck, ~var) );
        end
    end
    
    %% Split expr-unknown into expr-var and expr-expr
    for k=1:length(hsrf.eu)
        ck = hsrf.eu{k};
        [var, ID] = isa_variable_helper(ck, 'rhs');
        for u = yop.row_vec( unique(ID(var)) )
            hsrf.add_ev( yop.get_subrelation(ck, ID==u) );
        end
        
        if ~all(var)
            hsrf.add_ee( yop.get_subrelation(ck, ~var) );
        end
    end
    
    %% Reset unknown in order not to process any relation twice
    hsrf.clear_unknown();
end

end

function [vars, ID] = isa_variable_helper(relation, side)
% If one side is scalar, the other can be something different.
% The values are therefore scaled according to the size of the
% relation, which has a valid size.

[vars, ID] = isa_variable(relation.(side));

if isscalar(relation.(side))
    vars = vars & true(size(relation));
    ID = ID * ones(size(relation));
end
end

% function sr = get_subrelation(relation, idx)
% % Since indices might be scaled and variables can be scalars it is
% % necessary to test if it is possible to take the subindices of the
% % relations.
% 
% if all(idx==false)
%     sr = [];
%     return;
% end
% 
% if isscalar(relation.rhs) %&& ~all(idx==false)
%     rhs = relation.rhs;
% else
%     rhs = relation.rhs(idx);
% end
% 
% if isscalar(relation.lhs)  %&& ~all(idx==false)
%     lhs = relation.lhs;
% else
%     lhs = relation.lhs(idx);
% end
% 
% f  = get_constructor(relation);
% sr = f(lhs, rhs);
% 
% end