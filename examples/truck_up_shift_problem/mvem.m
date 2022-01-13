function [dx, y, c] = mvem(x, u, M_w, J_genset)
% MVEM2 - a diesel-electric powertrain model, with the engine modeled using
% typical engine efficiency characteristic as described in:
% "MODELLING FOR OPTIMAL CONTROL: A VALIDATED DIESEL-ELECTRIC POWERTRAIN MODEL"
% Martin Sivertsson and Lars Erisson
% Contact: marsi@isy.liu.se
%
%The  model has the states: x=[w_ice;p_im;p_em;w_tc]
% controls u=[u_f;u_wg;P_gen].
%MVEM2 takes x and u and param(the model parameters) and returns the
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
%
%     Modified 2019, Viktor Leek
%     Modified 2022, Viktor Leek
%
% -----------------------------------------------------------------------------
upshift_param;

w_ice = x(1);
p_im  = x(2);
p_em  = x(3);
w_tc  = x(4);

u_f   = u(1);
u_wg  = u(2);
P_gen = u(3);

%% Compressor
% Massflow
Pi_c = p_im / p_c_b;
w_tc_corr = w_tc / sqrt(T_c_b / T_amb);
w_tc_corr_norm = w_tc_corr / 15000;
Pi_c_max= ((w_tc_corr*R_c)^2 * Psi_max/(2*cp_a.*T_c_b) + 1)^(g_a/(g_a-1));
dot_m_c_corr_max = dot_m_c_corr_max(1)*(w_tc_corr_norm)^2 + ...
                   dot_m_c_corr_max(2)*(w_tc_corr_norm) + ...
                   dot_m_c_corr_max(3);
dot_m_c_corr = dot_m_c_corr_max * sqrt( 1 - (Pi_c / Pi_c_max).^2 );
dot_m_c = dot_m_c_corr * (p_c_b / p_amb) / sqrt(T_c_b /T_amb);

% Power
Phi = dot_m_c * R_a * T_c_b / (w_tc * 8 * R_c^3  * p_c_b);
dPhi = Phi - Phi_opt;
dN = w_tc_corr_norm - w_tc_corr_opt;
eta_c = eta_c_max - ( Q(1) * dPhi^2 + 2*dPhi*dN* Q(3) + Q(2) * dN^2 );
P_c = dot_m_c * cp_a * T_c_b * ...
     ( Pi_c^( (g_a - 1)/g_a) - 1)/eta_c;

%% Cylinder
% Airflow
eta_vol = c_eta_vol(1)*sqrt(p_im) + c_eta_vol(2)*sqrt(w_ice) + c_eta_vol(3);
dot_m_ci = eta_vol * p_im * w_ice * V_D / ( 4*pi * R_a * T_im);

% Fuelflow
dot_m_f = u_f * w_ice * n_cyl*(10^-6) /  ( 4*pi );
phi = (AFs * dot_m_f) / dot_m_ci;

% Torque
a = (  eta_ig_isl(2) - eta_ig_isl(1) ) / ( -eta_ig_isl(3)^2 );
b = -2 * a * eta_ig_isl(3);
eta_factor = a*( u_f / w_ice )^2 + b*u_f/w_ice + eta_ig_isl(1);
eta_ig = eta_factor * (1 - 1/( r_c^(g_cyl-1) ));
W_pump = V_D * (p_em - p_im);
W_ig = n_cyl * Hlhv * eta_ig * u_f * 10^-6;
W_fric = V_D * 10^5 * ( c_fr(1)*(w_ice)^2 + c_fr(2)*(w_ice) + c_fr(3) );
M_ice = ( W_ig - W_fric - W_pump) / (4*pi);
P_ice = M_ice * w_ice;

% Maximum fuel injection
u_f_max = 4*pi * 1e6 * dot_m_ci / ( w_ice * AFs * lambda_min * n_cyl );

%cylinder_TempOut
Pi_e = p_em / p_im;
q_in = dot_m_f * Hlhv / (dot_m_f + dot_m_ci);
T_eo = eta_sc * (Pi_e^(1-1/g_a)) * (r_c^(1-g_a)) * ...
    ( q_in / cp_a+T_im * r_c^(g_a-1) );
T_em = T_amb_r + (T_eo-T_amb_r) * ...
    exp( -h_tot*V_tot / ( (dot_m_f+dot_m_ci)*cp_e) );


%% turbine
%Massflow
Pi_t = p_es / p_em;
Psi_t = c_t(1) * sqrt( 1 - Pi_t^c_t(2) );
dot_m_t = p_em*Psi_t*A_t_eff/sqrt(T_em*R_e);

%Power
BSR = R_t * w_tc / sqrt( 2 * cp_e * T_em * (1 - Pi_t^(1 - 1/g_e)) );
eta_tm = eta_tm_max - c_m * (BSR - BSR_opt)^2;
P_t_eta_tm = dot_m_t * cp_e * T_em * eta_tm * (1 - Pi_t^((g_e-1)/g_e));

% wastegate_massflow
Psi_wg = c_wg(1) * sqrt(1 - Pi_t^c_wg(2));
dot_m_wg = p_em * Psi_wg * A_wg_eff * u_wg/sqrt(T_em * R_e);

%% Generator
a1 = gen2(1)*w_ice^2 + gen2(2)*w_ice + gen2(3);
a2 = gen2(4)*w_ice^2 + gen2(5)*w_ice + gen2(6);
f1 = a1*P_gen + gen2(7);
f2 = a2*P_gen + gen2(7);
g1 = ( 1 + tanh( 0.005 * P_gen ) )/2;
% P_mech = f2 + g1*( f1 - f2 );
P_mech = f2 + g1*( f1 - f2 ) - 3.60352247e2;
M_mech = P_mech/w_ice;

%% Dynamic eqs
dwice = (M_ice - M_mech - M_w)/J_genset;
dpim = R_a * T_im * ( dot_m_c - dot_m_ci ) / V_is;
dpem = R_e * T_em * (dot_m_ci + dot_m_f - dot_m_t - dot_m_wg) / V_em;
dwtc = (P_t_eta_tm - P_c) / (w_tc * J_tc);
dx = [dwice; dpim; dpem; dwtc];

%% Constraints
c = {};
c{end+1} = P_ice <= P_ice_max;
c{end+1} = P_ice <= cPice(1)*w_ice^2 + cPice(2)*w_ice + cPice(3);
c{end+1} = P_ice <= cPice(4)*w_ice^2 + cPice(5)*w_ice + cPice(6);
c{end+1} = phi <= 1/lambda_min;
c{end+1} = BSR_min <= BSR <= BSR_max;
c{end+1} = Pi_c <= c_mc_surge(1)*dot_m_c_corr + c_mc_surge(2);

%% Signals
y.compressor.speed = w_tc;
y.compressor.pressure_ratio = Pi_c;
y.compressor.efficiency = eta_c;
y.compressor.power = P_c;
y.compressor.surge_line = c_mc_surge(1) * dot_m_c_corr + c_mc_surge(2);

y.intake.pressure = p_im;
y.intake.temperature = T_im;

y.cylinder.volumetric_efficiency = eta_vol;
y.cylinder.air_massflow = dot_m_ci;
y.cylinder.fuel_injection = u_f;
y.cylinder.fuel_massflow = dot_m_f;
y.cylinder.fuel_to_air_ratio = phi;
y.cylinder.indicated_efficiency = eta_ig;
y.cylinder.indicated_torque = W_ig/(4*pi);
y.cylinder.pumping_torque = W_pump/(4*pi);
y.cylinder.friction_torque = W_fric/(4*pi);
y.cylinder.temperature_out = T_eo;
y.cylinder.fuel_max = u_f_max;
y.cylinder.lambda_min = 1.2;

y.phi_max = 1/1.2;
y.phi = phi;
y.u_f_max = u_f_max;

y.engine.speed = w_ice;
y.engine.efficiency = max(0, P_ice/(dot_m_f*Hlhv));
y.engine.torque = M_ice;
y.engine.power = P_ice;
y.engine.power_limit = [(cPice(1)*w_ice^2 + cPice(2)*w_ice + cPice(3)); ....
                       (cPice(4)*w_ice^2 + cPice(5)*w_ice + cPice(6))];

y.exhaust.pressure = p_em;
y.exhaust.temperature = T_em;

y.turbine.speed = w_tc;
y.turbine.pressure_ratio = Pi_t;
y.turbine.massflow = dot_m_t;
y.turbine.BSR = BSR;
y.turbine.BSR_max = BSR_max;
y.turbine.BSR_min = BSR_min;
y.turbine.efficiency = eta_tm;
y.turbine.power = P_t_eta_tm;

y.wastegate.control = u_wg;
y.wastegate.massflow = dot_m_wg;

y.turbocharger.speed = w_tc;

y.generator.power = P_gen;
y.generator.torque = P_gen/w_ice;
y.emachine.power = P_mech;
y.emachine.torque = M_mech;


end


