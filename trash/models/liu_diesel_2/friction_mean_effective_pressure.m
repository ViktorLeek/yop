function FMEP = friction_mean_effective_pressure(N_ice, u_f, param)
FMEP = param.c_fric(4)*N_ice.*u_f + param.c_fric(3)*u_f + param.c_fric(2)*N_ice + param.c_fric(1);
end

