function [W_thr, Psi, Pi] = W_throttle(pus, pds, T_c, u_thr, param)
Pi = pds./pus;
% Psi = Psi_thr(pus,pds,param);
Psi = Psi_ohata(Pi, param.gamma_air);
W_thr = pus./sqrt(param.R_a.*T_c) .* param.throttle.Cd .* param.throttle.A .* u_thr .* Psi;
end