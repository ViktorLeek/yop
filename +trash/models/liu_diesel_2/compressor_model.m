function [Wc,eta_c,T_out,Pc,h,Wch_SpL,PIch_SpL,WZSL_SpL,PIZSL_SpL] = compressor_model(p_in, p_out, T_in, omega_c, param)

%% Extract Model Parameters and preprocess the inputs
T_ref     	= param.T_ref;
p_ref     	= param.p_ref; 
Cp_air      = param.Cp_air; 
gamma_air   = param.gamma_air;

PiC         = p_out./p_in;
Nc_corr     = (omega_c.*(30./pi))./(sqrt(T_in./T_ref));  % corrected rad/s
Nc_corr_n   = Nc_corr./param.Nc_max_map;
% % Wc_corr     = zeros(length(PiC),1);

%% Mass flow model
% Compute base functions
Wch_SpL     = param.F_Wch(param.C_Wch,Nc_corr_n,param);
PIch_SpL    = param.F_PIch(param.C_PIch,Nc_corr_n,param);
WZSL_SpL    = param.F_WZS(param.C_Wzs,Nc_corr_n,param);
PIZSL_SpL   = param.F_PIZS(param.C_PIzs,Nc_corr_n,param);
CUR_SpL     = param.F_CUR(param.C_cur,Nc_corr_n);

% Calculate the ellipse value
Wc_corr 	= WZSL_SpL +(Wch_SpL-WZSL_SpL).*(1-((PiC-PIch_SpL)./(PIZSL_SpL-PIch_SpL)).^CUR_SpL).^(1./CUR_SpL);

% Uncorrect the massflow
Wc          = Wc_corr.*(p_in./p_ref)./(sqrt(T_in./T_ref));

%% Efficiency model 
% Compute base functions
B_model    	= param.F_b(param.C_b,Nc_corr_n,param);
A_model    	= -param.F_a(param.C_a,Nc_corr_n,param);
K_fric   	= param.F_kloss(param.C_loss,Nc_corr,Wc_corr,param);

% Compute actual Enthalphy
D_H_act   = K_fric.*(B_model+A_model.*Wc_corr);
% % D_H_act(D_H_act<1e-4) = 1e-4;

% Compute isentropic Enthalphy
Delta_h_is  = Cp_air.*T_ref.*(PiC.^((gamma_air - 1)./gamma_air) - 1);

% Compute efficiency with limitations
eta_c       = Delta_h_is./D_H_act;

%% Outlet temperature
T_out = T_in+(T_in./eta_c).*(PiC.^((gamma_air - 1)/gamma_air) - 1);

%% Consumed power
Pc  = Cp_air.*Wc.*(T_out-T_in);

%% Constraints
h = {...
%     PIch_SpL - PiC; ... % Onödigt om tryckkvoten måste vara större än 1.
    PiC - PIZSL_SpL; ... % Nödvändigt.
%     -D_H_act; ...
    -Delta_h_is; ... % Nödvändigt
%     eta_c-1; ...
%     -eta_c; ...
%     1-PiC; ...
    };
h = vertcat(h{:});


end
