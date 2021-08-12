function [T_t] = turbine_outlet_temperature(T_em, eta_tm, Pi_t, param)
    T_t = T_em - T_em .* eta_tm .* (1 - Pi_t .^ ((param.gamma_e - 1)/param.gamma_e));
end