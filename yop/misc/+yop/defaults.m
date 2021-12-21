classdef defaults
    properties (Constant) % OCP
        independent0_ub = 0;
        independent0_lb = 0;
        independentf_ub = inf;
        independentf_lb = 0;
        
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
        
        polynomial_degree = 5;
        collocation_points = 'legendre';
        control_invervals = 50;
        rk4_steps = 4;
    end
    
    properties (Constant) % IVP
        algebraic_guess = 1;
        ivp_solver = 'idas';
        ivp_sol_points = 100;
    end
    
    properties (Constant)
        independent_num = 0;
    end
end