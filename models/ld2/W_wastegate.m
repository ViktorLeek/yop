function [W_wg, Psi, Pi]= W_wastegate( pus, pds, T_em, u_wg, ice_param )
Pi = pds./pus;
% Psi = Psi_wg(pus,pds,ice_param);
Psi = Psi_ohata(Pi, ice_param.gamma_e);
W_wg = pus./sqrt(ice_param.R_e.*T_em) * ice_param.wastegate.A .* ice_param.wastegate.Cd .* u_wg .* Psi;
end