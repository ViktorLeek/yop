function [Wc,eta_c,T_out,Pc,h,Wch_SpL,PIch_SpL,WZSL_SpL,PIZSL_SpL] = compressor(p_in, p_out, T_in, omega_c, param)

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
Wch_SpL     = F_WCh(param.C_Wch, Nc_corr_n, param);
PIch_SpL    = F_PICh(param.C_PIch, Nc_corr_n, param);
WZSL_SpL    = F_WZS(param.C_Wzs, Nc_corr_n, param); % * FN
PIZSL_SpL   = F_PIZS(param.C_PIzs, Nc_corr_n, param);
CUR_SpL     = F_CUR(param.C_cur, Nc_corr_n);

%%
% Calculate the ellipse value
Wc_corr 	= WZSL_SpL +(Wch_SpL-WZSL_SpL).*(1-((PiC-PIch_SpL)./(PIZSL_SpL-PIch_SpL)).^CUR_SpL).^(1./CUR_SpL);

% Uncorrect the massflow
Wc          = Wc_corr.*(p_in./p_ref)./(sqrt(T_in./T_ref));

%% Efficiency model 
% Compute base functions
B_model    	= F_b_SAE(param.C_b,Nc_corr_n,param);
A_model    	= -F_a_SAE(param.C_a,Nc_corr_n,param);
K_fric   	= F_kloss(param.C_loss,Nc_corr,Wc_corr,param);

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

function [Wch_SpL] = F_WCh(c, Nc_corr_n, map)
Wch_SpL=map.Wc_max_map.*(c(1)+c(2).*atan(c(3).*Nc_corr_n-c(4)));
end

function [PIch_SpL] = F_PICh(c, Nc_corr_n, map)
PIch_SpL=map.PI_max_map.*(c(1)+c(2).*Nc_corr_n.^c(3));
end

function [WZSL_SpL] = F_WZS(c, Nc_corr_n, map)
WZSL_SpL=map.Wc_max_map.*(c(1).*Nc_corr_n.^c(2));
end

function [PIZSL_SpL] = F_PIZS(c, Nc_corr_n, map)
PIZSL_SpL   = 1+(map.PI_max_map-1).*(c(1).*Nc_corr_n.^c(2));
end

function [CUR_SpL] = F_CUR(c, Nc_corr_n)
CUR_SpL     = c(1)+c(2).*Nc_corr_n.^c(3);
end

function [a_p] = F_a_SAE(a,Nc_corr_n,map)
a_p = (map.Delta_h_max./map.Wc_max_map).*(a(1).*Nc_corr_n./((1+a(2).*Nc_corr_n.^2).^a(3)));
end

function [b_p] = F_b_SAE(c,Nc_corr_n,map)
b_p = map.Delta_h_max.*(c(1).*Nc_corr_n.^2+c(2).*Nc_corr_n.^3);
end

function [k_loss] = F_kloss(c,Nc_corr,Wc,map)
k_loss =1+c.*(map.rho1.*map.D_2.^3.*pi./60).*(Nc_corr./Wc);
end