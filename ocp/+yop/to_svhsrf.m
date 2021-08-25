function svhsrf = to_svhsrf(hsrf, svhsrf)
% converts the relations to single variable homogeneous single
% relation form. This means that relations that are of the form
% var (R) expr, or expr (R) var, is separated so that the
% variable side only contains elements from a single
%'yop.ast_variable'.

if nargin == 1
    svhsrf = yop.svhsrf_data();
end

for k=1:length(hsrf.ve)
    rk = hsrf.ve{k};
    constructor = yop.get_constructor(rk);
    
    % In order to split the variables it is first necessary to
    % deterimine the variables that reaches the definition of
    % lhs
    rv = reaching_variables(rk.lhs);
    
    % Based on the reaching variables the expression is divided
    % into subexpressions.
    for n=1:length(rv)
        if ~isempty(rv(n).reaching)
            idx = rv(n).index;
            if isscalar(rk.rhs)
                node = yop.ast_le(rk.lhs(idx), rk.rhs);
            else
                node = yop.ast_le(rk.lhs(idx), rk.rhs(idx));
            end
            svhsrf.add_ve(node);
        end
    end
    assert(~isempty(svhsrf.ve), ...
        ['[Yop] Unexpected error: Some variables are', ...
        'expected to reach the expression. Check if the', ...
        'proper transforms has been applied first.'])
end

for k=1:length(hsrf.ev)
    rk = hsrf.ev{k};
    constructor = yop.get_constructor(rk);
    
    % In order to split the variables it is first necessary to
    % deterimine the variables that reaches the definition of
    % lhs
    rv = reaching_variables(rk.rhs);
    
    % Based on the reaching variables the expression is divided
    % into subexpressions.
    for n=1:length(rv)
        if ~isempty(rv(n).reaching)
            idx = rv(n).index;
            if isscalar(rk.lhs)
                node = yop.ast_le(rk.lhs, rk.rhs(idx));
            else
                node = yop.ast_le(rk.lhs(idx), rk.rhs(idx));
            end
            svhsrf.add_ev(node);
        end
    end
    assert(~isempty(svhsrf.ev), ...
        ['[Yop] Unexpected error: Some variables are', ...
        'expected to reach the expression. Check if the', ...
        'proper transforms has been applied first.'])
end
end