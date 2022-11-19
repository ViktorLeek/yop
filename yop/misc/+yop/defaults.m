classdef defaults
    properties (Constant) % OCP
        independent0_ub = inf;
        independent0_lb = -inf;
        independentf_ub = inf;
        independentf_lb = -inf;
        
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
        
        state_der_ub = inf;
        state_der_lb = -inf;
        
        solver = 'ipopt'
        state_degree = 5;
        state_points = {'legendre'};
        control_degree = 0;
        control_points = {'radau'};
        control_invervals = 50;
        continuity = false;
        rk4_steps = 4;
    end
    
    properties (Constant) % IVP
        algebraic_guess = 1;
        ivp_solver = 'ode15s';
        ivp_sol_points = 100;
    end
    
    properties (Constant)
        independent_num = 0;
    end
end