function [ M_ice, aux ] = engine_torque(N_ice, p_im, p_em, u_f, phi, param)

% Ne_rps = N_ice/60;    % -> rps
% w_ice  = N_ice*pi/30; % -> rad/s

%% Indicated work
[eta_ig, eta_otto, gamma, eta_cal] = indicated_efficiency(phi, N_ice, u_f, param);
W_ig = param.q_HV .* eta_ig .* u_f * 1e-6 * param.n_cyl;
IMEP = W_ig./param.V_D;

%% Friction work
FMEP = friction_mean_effective_pressure(N_ice, u_f, param);
W_fric = param.V_D.*FMEP;

%% Pump work
PMEP = pump_mean_effective_pressure(p_em, p_im, param);
W_pump = param.V_D .* PMEP;

%% Engine Toruqe
M_ice = (W_ig - W_pump - W_fric)./(4*pi);
BMEP = M_ice*2*pi.*param.n_r./param.V_D;

%% Submodel signals
aux = struct;
aux.eta_ig = eta_ig;
aux.eta_otto = eta_otto;
aux.gamma = gamma;
aux.eta_cal = eta_cal;
aux.W_pump = W_pump;
aux.W_ig = W_ig;
aux.W_fric = W_fric;
aux.M_pump = W_pump/4/pi;
aux.M_ig = W_ig/4/pi;
aux.M_fric = W_fric/4/pi;
aux.BMEP = BMEP;
aux.IMEP = IMEP;
aux.PMEP = PMEP;
aux.FMEP = FMEP;

end

