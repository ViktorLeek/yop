function BSR = turbine_BSR(Pi_t, w_tc, T_em, param)

BSR = param.turb.R_t*w_tc./sqrt( 2*param.c_pe*T_em.*(1 - Pi_t.^(1-1/param.gamma_e)) );


end