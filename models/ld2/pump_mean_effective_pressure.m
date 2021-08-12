function PMEP = pump_mean_effective_pressure(p_em, p_im, param)
PMEP = param.c_pmep(1) + param.c_pmep(2)*(p_em - p_im);
end