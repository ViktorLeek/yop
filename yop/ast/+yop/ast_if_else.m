classdef ast_if_else < yop.ast_expression
    properties
        m_cond
        m_true
        m_false
        m_short % short circuit option
    end
    methods
        function obj = ast_if_else(cond, True, False, short)
            if nargin == 3
                short = true;
            end
            
            sz = size(True);
            val = if_else(value(cond), value(True), value(False), short);
            num = nan(sz);
            
            if ~isinf(get_t0(True)) && ~isinf(get_t0(False))
                t0 = max(get_t0(True), get_t0(False));
            elseif ~isinf(get_t0(True))
                t0 = get_t0(True);
            elseif ~isinf(get_t0(False))
                t0 = get_t0(False);
            else
                t0 = yop.initial_timepoint*ones(sz);
            end
            
            if ~isinf(get_tf(True)) && ~isinf(get_tf(False))
                tf = min(get_tf(True), get_tf(False));
            elseif ~isinf(get_tf(True))
                tf = get_tf(True);
            elseif ~isinf(get_tf(False))
                tf = get_tf(False);
            else
                tf = yop.final_timepoint*ones(sz);
            end
            
            if isa_reducible(True) && isa_reducible(False)
                red = true(sz);
            else
                red = false(sz);
            end
            
            obj@yop.ast_expression( ...
                val      , ... value
                num      , ... numval
                t0       , ... t0
                tf       , ... tf
                false(sz), ... der
                red      , ... reducible
                zeros(sz), ... type
                zeros(sz) ... typeid
                );
            
            obj.m_cond = cond;
            obj.m_true = True;
            obj.m_false = False;
            obj.m_short = short;
        end
    end
end