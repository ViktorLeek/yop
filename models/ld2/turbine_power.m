function P_turbine = turbine_power(W_turbine, T_em, PiT, eta_tm, param)

P_turbine = W_turbine*param.c_pe.*T_em.*eta_tm.*(1-PiT.^(1-1/param.gamma_e));

end