function [eta_t, BSR_opt] = turbine_efficiency(BSR, N_tc_corr, ice_param)

N_tc_norm = N_tc_corr.*1e-4;

eta_t_max   = ice_param.turb.c_eta0 + ice_param.turb.c_eta1.*N_tc_norm + ice_param.turb.c_eta2.*N_tc_norm.^2;

k_eta       = ice_param.turb.c_max0 + ice_param.turb.c_max1.*N_tc_norm;

BSR_opt     = ice_param.turb.c_BSR0 + ice_param.turb.c_BSR1.*N_tc_norm + ice_param.turb.c_BSR2.*N_tc_norm.^2;

eta_t       = eta_t_max + k_eta.*(BSR-BSR_opt).^2;

end

