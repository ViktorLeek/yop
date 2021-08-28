function [W_t_corr, k0, k1, k2] = turbine_massflow(PiT, N_t_corr, param)


N_t_norm = N_t_corr./1000;

k00 = param.turb.k0(1);
k02 = param.turb.k0(2);

k10 = param.turb.k1(1);
k11 = param.turb.k1(2);

k20 = param.turb.k2(1);
k21 = param.turb.k2(2);
k22 = param.turb.k2(3);

k0 = k02*N_t_norm.^2                + k00;
k1 =                   k11*N_t_norm + k10;
k2 = k22*N_t_norm.^2 + k21*N_t_norm + k20;

W_t_corr = k0.*(1 - PiT.^k1).^k2;

end