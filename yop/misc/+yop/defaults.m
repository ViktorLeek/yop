classdef defaults
    properties (Constant)
        
        independent_ub0 = 0;
        independent_lb0 = 0;
        independent_ub  = inf;  % semantics. t = [lb, ub], fixed time?
        independent_lb  = 0;
        independent_ubf = inf;
        independent_lbf = 0;
        
        state_ub0 = inf;
        state_ubf = inf;
        state_ub  = inf;
        state_lb0 = -inf;
        state_lbf = -inf;
        state_lb  = -inf;
        
        algebraic_ub = inf;
        algebraic_lb = -inf;
        
        control_ub0 = inf;
        control_ubf = inf;
        control_ub  = inf;
        control_lb0 = -inf;
        control_lbf = -inf;
        control_lb  = -inf;
        
        parameter_ub = inf;
        parameter_lb = -inf;
        
    end
end