function [dx,c]=MVEMo(w_ice,p_im,p_em,w_tc,u_f,u_wg,P_gen,param)
% MVEMo - a diesel-electric powertrain model, with the engine modeled using
% real engine efficiency characteristic as described in:
% "MODELLING FOR OPTIMAL CONTROL: A VALIDATED DIESEL-ELECTRIC POWERTRAIN MODEL"
% Martin Sivertsson and Lars Erisson
% Contact: marsi@isy.liu.se
% 
%The  model has the states: x=[w_ice;p_im;p_em;w_tc]
% controls u=[u_f;u_wg;P_gen]. 
%MVEMo takes x and u and param(the model parameters) and returns the
%state derivatives: dx=[dwice;dpim;dpem;dwtc] and some additional variables
%in the struct c.
%
%-----------------------------------------------------------------------------
%     Copyright 2014, Martin Sivertsson, Lars Eriksson
% 
%     This package is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as
%     published by the Free Software Foundation, version 3 of the License.
% 
%     This package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------

W_ICE=param.state_norm(1)*w_ice;
P_IM=param.state_norm(2)*p_im;
P_EM=param.state_norm(3)*p_em;
W_TC=param.state_norm(4)*w_tc;

U_F=param.control_norm(1)*u_f;
P_GEN=param.control_norm(3)*P_gen;







%Compressor
% Massflow
Pi_c=P_IM./param.p_c_b;
w_tc_corr=W_TC./sqrt(param.T_c_b./param.T_amb);
w_tc_corr_norm=w_tc_corr/15000;
Pi_c_max=(((w_tc_corr.*param.R_c).^2)*param.Psi_max./(2*param.cp_a.*param.T_c_b)+1).^(param.gamma_a/(param.gamma_a-1));
dot_m_c_corr_max=param.dot_m_c_corr_max(1)*(w_tc_corr_norm).^2+param.dot_m_c_corr_max(2).*(w_tc_corr_norm)+param.dot_m_c_corr_max(3);
dot_m_c_corr=dot_m_c_corr_max.*(sqrt(1-(Pi_c./Pi_c_max).^2));
dot_m_c=dot_m_c_corr.*(param.p_c_b./param.p_amb)./sqrt(param.T_c_b./param.T_amb);
% Power
Phi=dot_m_c.*param.R_a.*param.T_c_b./(W_TC.*8.*param.R_c.^3.*param.p_c_b);
dPhi=Phi-param.Phi_opt;
dN=(w_tc_corr_norm)-param.w_tc_corr_opt;
eta_c=param.eta_c_max-(param.Q(1)*dPhi.^2+2*dPhi.*dN* param.Q(3)+param.Q(2)*dN.^2);
P_c=dot_m_c.*param.cp_a.*param.T_c_b.*(Pi_c.^((param.gamma_a-1)/param.gamma_a)-1)./eta_c;

%cylinder_Airflow
eta_vol=param.c_eta_vol(1)*sqrt(P_IM)+param.c_eta_vol(2)*sqrt(W_ICE)+param.c_eta_vol(3);
dot_m_ci=eta_vol.*P_IM.*W_ICE.*param.V_D./(4*pi*param.R_a*param.T_im);

% cylinder_Fuelflow
dot_m_f=U_F.*W_ICE*param.n_cyl*(10^-6)/(4*pi);
%lambda=dot_m_ci./(param.AFs*dot_m_f); %Becomes infinite for zero massflow
phi=(param.AFs*dot_m_f)./dot_m_ci;

%cylinder_Torque
gf=(1+tanh(0.1*(W_ICE-1500*pi/30)))/2;
Mf1=param.eta_ig_w2(1)*W_ICE.^2+param.eta_ig_w2(2)*W_ICE;
Mf2=param.eta_ig_w2(3)*W_ICE.^2+param.eta_ig_w2(4)*W_ICE+param.eta_ig_w2(5);
eta_factor=Mf1+gf.*(Mf2-Mf1);
eta_ig=eta_factor*(1-1/(param.r_c^(param.gamma_cyl-1)));
M_pump=param.V_D*(P_EM-P_IM);
M_ig=param.n_cyl*param.Hlhv*eta_ig.*U_F*(10^-6);
M_fric=param.V_D*(10^5)*(param.c_fr(1)*(W_ICE).^2+param.c_fr(2)*(W_ICE)+param.c_fr(3));
M_ice=(M_ig-M_fric-M_pump)/(4*pi);
P_ice=M_ice.*W_ICE;


%cylinder_TempOut
Pi_e=P_EM./P_IM;
q_in=dot_m_f.*param.Hlhv./(dot_m_f+dot_m_ci);
T_eo=param.eta_sc*(Pi_e.^(1-1/param.gamma_a)).*(param.r_c.^(1-param.gamma_a)).*(q_in./param.cp_a+param.T_im*param.r_c^(param.gamma_a-1));
T_em=param.T_amb_r+(T_eo-param.T_amb_r).*exp(-param.h_tot*param.V_tot./((dot_m_f+dot_m_ci)*param.cp_e));


% turbine
%Massflow
Pi_t=param.p_es./P_EM;
Psi_t=param.c_t(1)*sqrt(1-power(Pi_t,param.c_t(2)));
dot_m_t=P_EM.*Psi_t*param.A_t_eff./sqrt(T_em*param.R_e);
%Power
BSR=param.R_t*W_TC./(sqrt(2*param.cp_e*T_em.*(1-Pi_t.^(1-1/param.gamma_e))));
eta_tm=param.eta_tm_max-param.c_m.*(BSR-param.BSR_opt).^2;
P_t_eta_tm=dot_m_t.*param.cp_e.*T_em.*eta_tm.*(1-Pi_t.^((param.gamma_e-1)/param.gamma_e));


% wastegate_massflow
Psi_wg=param.c_wg(1)*sqrt(1-power(Pi_t,param.c_wg(2)));
dot_m_wg=P_EM.*Psi_wg.*param.A_wg_eff.*u_wg./sqrt(T_em*param.R_e);

%generator
a1=param.gen2(1).*W_ICE.^2+param.gen2(2).*W_ICE+param.gen2(3);
a2=param.gen2(4).*W_ICE.^2+param.gen2(5).*W_ICE+param.gen2(6);
f1=a1.*P_GEN+param.gen2(7);
f2=a2.*P_GEN+param.gen2(7);
g1=(1+tanh(0.005*P_GEN))./2;
P_mech=f2+g1.*(f1-f2);

%Dynamic equations
dwice=(P_ice-P_mech)./W_ICE./param.J_genset/param.state_norm(1);
dpim=param.R_a*param.T_im*(dot_m_c-dot_m_ci)/param.V_is/param.state_norm(2);
dpem=param.R_e*T_em.*(dot_m_ci+dot_m_f-dot_m_t-dot_m_wg)/param.V_em/param.state_norm(3);
dwtc=(((P_t_eta_tm-P_c)./(W_TC*param.J_tc)))/param.state_norm(4);

dx = [dwice; dpim; dpem; dwtc];

c=struct('W_ICE',W_ICE,'P_IM',P_IM,'P_EM',P_EM,'W_TC',W_TC,...
    'U_F',U_F,'U_WG',u_wg,'P_GEN',P_GEN,'dot_m_f',dot_m_f,...
    'Pi_c',Pi_c,'P_ice',P_ice,'BSR',BSR,'phi',phi,'dot_m_c_corr',dot_m_c_corr);

end


