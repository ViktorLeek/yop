%% -------- Transient optimization of Diesel-electric Powertrain ----------
%
%  Copyright 2019, Viktor Leek
%                  viktor.leek@liu.se
%
% -------------------------------------------------------------------------
% Problem taken from:
% AN OPTIMAL CONTROL BENCHMARK: TRANSIENT OPTIMIZATIONOF A DIESEL-ELECTRIC
%   POWERTRAIN
%
% Formulated by: Martin Sivertsson and Lars Eriksson
% Publication:
%        http://www.fs.isy.liu.se/Publications/Articles/SIMS_14_MS_LE_2.pdf
%
% Software available at:
%        http://www.fs.isy.liu.se/Software/LiU-D-El_and_Benchmark/
%
% -------------------------------------------------------------------------
%% System of interest
genset = YopSystem('states', 5, 'controls', 3, 'model', @gensetModel);

time         = genset.t;
engine       = genset.y.engine;
compressor   = genset.y.compressor;
intake       = genset.y.intake;
cylinder     = genset.y.cylinder;
exhaust      = genset.y.exhaust;
turbine      = genset.y.turbine;
wastegate    = genset.y.wastegate;
turbocharger = genset.y.turbocharger;
generator    = genset.y.generator;

%% Initial Guess
% An engine speed controller
desiredEngineSpeed = 1500; % rpm
controller = YopSystem('states', 1, 'externals', 3, 'parameters', 2, 'model', @speedController);

% A ramp for controlling the generator power output
rampStartingPoint = 1;
rampFinishPoint = 4;
preRampValue = 0;
postRampValue = 120e3;
rampParameters = [rampStartingPoint; rampFinishPoint; preRampValue; postRampValue];

demand = YopSystem('model', @(t) powerDemand(t, rampParameters) );

% Connect the systems
closed = 0;
c1 = YopConnection(controller.y.engineSpeedInput, engine.speed);
c2 = YopConnection(controller.y.desiredEngineSpeedInput, rpm2rad(desiredEngineSpeed));
c3 = YopConnection(controller.y.fuelLimiterInput, cylinder.fuelLimiter);
c4 = YopConnection(controller.y.controlSignal, cylinder.fuelInjection);
c5 = YopConnection(wastegate.control, closed);
c6 = YopConnection(generator.power, demand.y.power);

% Create a simulator
simulator = YopSimulator(...
    'systems', [genset; controller; demand], ...
    'connections', [c1; c2; c3; c4; c5; c6] ...
    );

% Simulate the system
res = simulator.simulate(...
    'grid', linspace(0, 7, 1000), ...
    'printStats', true, ...
    'initialValue', engine.speed, rpm2rad(800), ...
    'initialValue', intake.pressure, 1.0143e+05, ...
    'initialValue', exhaust.pressure, 1.0975e+05, ...
    'initialValue', turbocharger.speed, 2.0502e+03, ...
    'initialValue', generator.energy, 0, ...
    'initialValue', controller.y.integralState, 0, ...
    'initialValue', controller.y.proportionalGain, 2, ...
    'initialValue', controller.y.integralGain, 1 ...
    );

%% Plot initial guess
% States
figure(1)
ax1 = subplot(511); hold on; grid on
res.plot(time, rad2rpm(engine.speed))
title('Engine speed [rpm]')

ax2 = subplot(512); hold on; grid on
res.plot(time, intake.pressure*1e-5)
title('Intake manifold pressure [bar]')

ax3 = subplot(513); hold on; grid on
res.plot(time, exhaust.pressure*1e-5)
title('Exhaust manifold pressure [bar]')

ax4 = subplot(514); hold on; grid on
res.plot(time, rad2rpm(turbocharger.speed)*1e-3)
title('Turbo speed [krpm]')

ax5 = subplot(515); hold on; grid on
res.plot(time, generator.energy*1e-3)
title('Generator energy output [kJ]')

% Controls
figure(2)
ax6 = subplot(311); hold on; grid on
res.plot(time, cylinder.fuelInjection);
res.plot(time, cylinder.fuelLimiter, 'r')
title('Fuel injection [mg/cycle/cylinder]')
legend('u_{fuel}', 'smoke limiter', 'Location', 'South')

ax7 = subplot(312); hold on; grid on
res.plot(time, wastegate.control);
title('Wastegate range[0,1]')

ax8 = subplot(313); hold on; grid on
res.plot(time, generator.power*1e-3);
title('Generator power [kW]')

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')


%% Optimal control


ocp = YopOcp();
ocp.min({ timeIntegral( cylinder.fuelMassflow ) })
ocp.st(...
    'systems', genset, ...
    { 1.1   '<=' t_f( time ) '<=' 1.4   }, ...
    ... Initial conditions
    { rpm2rad(800) '==' t_0( engine.speed      )  }, ...
    { 1.0143e+05   '==' t_0( intake.pressure   )  }, ...
    { 1.0975e+05   '==' t_0( exhaust.pressure  )  }, ...
    { 2.0502e+03   '==' t_0( turbocharger.speed)  }, ...
    {     0        '==' t_0( generator.energy  )  }, ...
    {  closed      '==' t_0( wastegate.control )  }, ...
    ... Terminal conditions
    { 100e3 '==' t_f( generator.power  ) }, ...
    { 100e3 '==' t_f( generator.energy ) }, ...
    {  0    '==' t_f( genset.ode(1:4)  ) }, ... % Stationarity
    ... Box constraints
    { rpm2rad(800) '<=' engine.speed           '<=' rpm2rad(2500) }, ...
    {  8.0889e+04  '<=' intake.pressure        '<=' 350000        }, ...
    {  9.1000e+04  '<=' exhaust.pressure       '<=' 400000        }, ...
    {     500      '<=' turbocharger.speed     '<=' 15000         }, ...
    {      0       '<=' generator.energy       '<=' 3000000       }, ...
    {      0       '<=' cylinder.fuelInjection '<=' 150           }, ...
    {      0       '<=' wastegate.control      '<=' 1             }, ...
    {      0       '<=' generator.power        '<=' 100e3         }, ...
    ... Path constraints
    { turbine.BSRMin '<=' turbine.BSR '<=' turbine.BSRMax         }, ...
    {  0 '>=' (engine.power - engine.powerLimit(1))               }, ...
    {  0 '>=' (engine.power - engine.powerLimit(2))               }, ...
    {  0 '>=' (cylinder.fuelToAirRatio - 1/cylinder.lambdaMin)    }, ...
    {  0 '>=' (compressor.pressureRatio - compressor.surgeline)   } ...
    );

% Scale problem
ocp.scale('objective', 1e3);

ocp.scale('variable', engine.speed,      'weight', rpm2rad(1e-3))
ocp.scale('variable', intake.pressure,   'weight', 1e-5)
ocp.scale('variable', exhaust.pressure,  'weight', 1e-5)
ocp.scale('variable', turbocharger.speed,'weight', 1e-3)
ocp.scale('variable', generator.energy,  'weight', 1e-5)

ocp.scale('variable', cylinder.fuelInjection, 'weight', 1e-2)
ocp.scale('variable', wastegate.control, 'weight', 1e-5)
ocp.scale('variable', generator.power,   'weight', 1e-5)

% Solve
sol = ocp.solve( ...
    'initialGuess', res, ...
    'controlIntervals', 75, ...
    'collocationPoints', 'radau', ...
    'polynomialDegree', 3, ...
    'ipopt', struct('acceptable_tol', 1e-7) ...
    );

%% Plot optimal control
% States
figure(1)
ax1 = subplot(511); hold on; grid minor
sol.plot(time, rad2rpm(engine.speed))
title('Engine speed [rpm]')

ax2 = subplot(512); hold on; grid minor
sol.plot(time, intake.pressure*1e-5)
title('Intake manifold pressure [bar]')

ax3 = subplot(513); hold on; grid minor
sol.plot(time, exhaust.pressure*1e-5)
title('Exhaust manifold pressure [bar]')

ax4 = subplot(514); hold on; grid minor
sol.plot(time, rad2rpm(turbocharger.speed)*1e-3)
title('Turbo speed [krpm]')

ax5 = subplot(515); hold on; grid minor
sol.plot(time, generator.energy*1e-3)
title('Generator energy output [kJ]')

% Controls
figure(2)
ax6 = subplot(311); hold on; grid minor
sol.stairs(time, cylinder.fuelInjection);
sol.plot(time, cylinder.fuelLimiter, 'r')
title('Fuel injection [mg/cycle/cylinder]')
legend('u_{fuel}', 'smoke limiter', 'Location', 'South')



ax7 = subplot(312); hold on; grid minor
sol.stairs(time, wastegate.control);
title('Wastegate range[0,1]')

ax8 = subplot(313); hold on; grid minor
sol.stairs(time, generator.power*1e-3);
title('Generator power [kW]')

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')

%% Functions
function out = rpm2rad(in)
out = in*pi/30;
end

function out = rad2rpm(in)
out = in*30/pi;
end

function y = powerDemand(t, p)
%% ------- Generator power demand -----------
%
%  Author: Viktor Leek, viktor.leek@liu.se
%  Copyright 2019
%
%--------------------------------------------
import casadi.*

rampStartingPoint = p(1);
rampFinishPoint = p(2);
preRampValue = p(3);
postRampValue = p(4);

tVec     = [   0,   rampStartingPoint, rampFinishPoint, 1000];
valueVec = [preRampValue, preRampValue, postRampValue, postRampValue];
lookupTable = casadi.interpolant('power', 'linear', {tVec}, valueVec);

y.power = lookupTable(t);
y.rampStartingPoint = rampStartingPoint;
y.rampFinishPoint = rampFinishPoint;
y.preRampValue = preRampValue;
y.postRampValue = postRampValue;
end


function [stateDerivative, y] = speedController(time, integralState, external, parameters)
%% A genset engine speed controller ---------
%  PI-controller with wind-up (!)
%
%  Author: Viktor Leek, viktor.leek@liu.se
%  Copyright 2019
%
%--------------------------------------------
engineSpeed = external(1); % Engine speed
desiredEngineSpeed = external(2); % Desired engine speed
fuelLimiter = external(3); % Smoke limiter

% Controller parameters
kp = parameters(1);
ki = parameters(2);

error = desiredEngineSpeed - engineSpeed;

% Feedback term
proportionalControl = kp*error;
integralControl = ki*integralState;
control = proportionalControl + integralControl;

% derivative of the integral
stateDerivative = error;

% Saturate contol signal
saturatedControl = if_else(control > fuelLimiter, fuelLimiter, control);

% Absolute limit
controlUpperBounded = if_else(saturatedControl < 0, 0, saturatedControl);
controlUpperLowerBounded = if_else(controlUpperBounded > 150, ...
                                    150, ...
                                    controlUpperBounded ...
                                   );

% Store output as a struct:
y.engineSpeedInput = external(1);
y.desiredEngineSpeedInput = external(2);
y.fuelLimiterInput = external(3);
y.proportionalGain = parameters(1);
y.integralGain = parameters(2);
y.integralState = integralState;
y.controlSignal = controlUpperLowerBounded; % Engine control signal
y.unsaturatedControl = control; % unconstrained control signal
y.error = error; % Error
y.integralControl = integralControl; % integral part
y.proportionalControl = proportionalControl; % proportional part
y.fuelLimiter = fuelLimiter; % Smoke limiter
y.externalInput = external;
end

function [dX, y] = gensetModel(time, state, control)
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
%
% -----------------------------------------------------------------------------
w_ice = state(1);
p_im = state(2);
p_em = state(3);
w_tc = state(4);
E_gen = state(5);

u_f = control(1);
u_wg = control(2);
P_gen = control(3);

param = gensetParameters;

%% Compressor
% Massflow
Pi_c = p_im / param.p_c_b;

w_tc_corr = w_tc / sqrt(param.T_c_b / param.T_amb);

w_tc_corr_norm = w_tc_corr / 15000;

Pi_c_max=( ...
           (w_tc_corr*param.R_c)^2 * param.Psi_max/(2*param.cp_a.*param.T_c_b) + 1 ...
         )^( param.gamma_a / (param.gamma_a-1) );

dot_m_c_corr_max = param.dot_m_c_corr_max(1)*(w_tc_corr_norm)^2 + ...
                   param.dot_m_c_corr_max(2)*(w_tc_corr_norm) + ...
                   param.dot_m_c_corr_max(3);

dot_m_c_corr = dot_m_c_corr_max * sqrt( 1 - (Pi_c / Pi_c_max).^2 );

dot_m_c = dot_m_c_corr * (param.p_c_b / param.p_amb) / ...
                     sqrt(param.T_c_b /param.T_amb);

% Power
Phi = dot_m_c * param.R_a * param.T_c_b / ...
     (w_tc * 8 * param.R_c^3  * param.p_c_b);

dPhi = Phi - param.Phi_opt;

dN = w_tc_corr_norm - param.w_tc_corr_opt;

eta_c = param.eta_c_max - ...
    ( param.Q(1) * dPhi.^2 + 2*dPhi.*dN* param.Q(3) + ...
      param.Q(2) * dN.^2 );

P_c = dot_m_c * param.cp_a * param.T_c_b * ...
     ( Pi_c^( (param.gamma_a - 1)/param.gamma_a) - 1)/eta_c;

%% Cylinder
% Airflow
eta_vol = param.c_eta_vol(1) * sqrt(p_im) + ....
          param.c_eta_vol(2) * sqrt(w_ice) + ...
          param.c_eta_vol(3);

dot_m_ci = eta_vol * p_im * w_ice * param.V_D / ...
          ( 4*pi * param.R_a * param.T_im);

% Fuelflow
dot_m_f = u_f * w_ice * param.n_cyl*(10^-6) / ...
         ( 4*pi );

phi = (param.AFs * dot_m_f) / dot_m_ci;

% Torque
a = (  param.eta_ig_isl(2) - param.eta_ig_isl(1) ) / ...
    ( -param.eta_ig_isl(3)^2 );

b = -2 * a * param.eta_ig_isl(3);

eta_factor = a*( u_f / w_ice )^2 + b*u_f/w_ice + param.eta_ig_isl(1);

eta_ig = eta_factor * (1 - 1/( param.r_c^(param.gamma_cyl-1) ));

W_pump = param.V_D * (p_em - p_im);

W_ig = param.n_cyl * param.Hlhv * eta_ig * u_f * 10^-6;

W_fric = param.V_D * 10^5 * ...
    ( param.c_fr(1)*(w_ice).^2 + ...
      param.c_fr(2)*(w_ice) + ...
      param.c_fr(3) );

M_ice = ( W_ig - W_fric - W_pump) / (4*pi);

P_ice = M_ice * w_ice;

% Maximum fuel injection
u_f_max = 4*pi * 1e6 * dot_m_ci / ...
    ( w_ice * param.AFs * param.lambda_min * param.n_cyl );

%cylinder_TempOut
Pi_e = p_em / p_im;

q_in = dot_m_f * param.Hlhv / (dot_m_f + dot_m_ci);

T_eo = param.eta_sc * (Pi_e^( 1 - 1/param.gamma_a )) * ...
      (param.r_c.^(1-param.gamma_a)) * ...
      ( q_in / param.cp_a+param.T_im * param.r_c^(param.gamma_a-1) );

T_em = param.T_amb_r + (T_eo-param.T_amb_r) * ...
    exp( -param.h_tot*param.V_tot / ( (dot_m_f+dot_m_ci)*param.cp_e) );


%% turbine
%Massflow
Pi_t = param.p_es / p_em;

Psi_t = param.c_t(1) * sqrt( 1 - Pi_t^param.c_t(2) );

dot_m_t=p_em.*Psi_t*param.A_t_eff./sqrt(T_em*param.R_e);

%Power
BSR = param.R_t * w_tc / ...
      sqrt( 2 * param.cp_e * T_em * (1 - Pi_t^(1 - 1/param.gamma_e)) );

eta_tm = param.eta_tm_max - param.c_m * (BSR - param.BSR_opt)^2;

P_t_eta_tm = dot_m_t * param.cp_e * T_em * eta_tm * ...
            (1 - Pi_t^( (param.gamma_e - 1) / param.gamma_e) );

% wastegate_massflow
Psi_wg = param.c_wg(1) * sqrt(1 - Pi_t^param.c_wg(2));

dot_m_wg = p_em * Psi_wg * param.A_wg_eff * u_wg/sqrt(T_em * param.R_e);

%% Generator
a1 = param.gen2(1)*w_ice^2 + param.gen2(2)*w_ice + param.gen2(3);
a2 = param.gen2(4)*w_ice^2 + param.gen2(5)*w_ice + param.gen2(6);

f1 = a1*P_gen + param.gen2(7);
f2 = a2*P_gen + param.gen2(7);

g1 = ( 1 + tanh( 0.005 * P_gen ) )/2;

P_mech = f2 + g1*( f1 - f2 );

%% Dynamic equations
dwice = (P_ice - P_mech) / w_ice / param.J_genset;

dpim = param.R_a * param.T_im * ( dot_m_c - dot_m_ci ) / param.V_is;

dpem = param.R_e * T_em * (dot_m_ci + dot_m_f - dot_m_t - dot_m_wg) / param.V_em;

dwtc = (P_t_eta_tm - P_c) / (w_tc * param.J_tc);

dE_gen = P_gen;

dX = [dwice; dpim; dpem; dwtc; dE_gen];

%% Signals
y.compressor.speed = w_tc;
y.compressor.pressureRatio = Pi_c;
y.compressor.efficiency = eta_c;
y.compressor.power = P_c;
y.compressor.surgeline = param.c_mc_surge(1) * dot_m_c_corr + param.c_mc_surge(2);

y.intake.pressure = p_im;
y.intake.temperature = param.T_im;

y.cylinder.volumetricEfficiency = eta_vol;
y.cylinder.airMassflow = dot_m_ci;
y.cylinder.fuelInjection = u_f;
y.cylinder.fuelMassflow = dot_m_f;
y.cylinder.fuelToAirRatio = phi;
y.cylinder.indicatedEfficiency = eta_ig;
y.cylinder.indicatedTorque = W_ig/(4*pi);
y.cylinder.pumpingTorque = W_pump/(4*pi);
y.cylinder.frictionTorque = W_fric/(4*pi);
y.cylinder.temperatureOut = T_eo;
y.cylinder.fuelLimiter = u_f_max;
y.cylinder.lambdaMin = 1.2;

y.engine.speed = w_ice;
y.engine.efficiency = if_else(u_f <= 0, 0, P_ice/(dot_m_f*param.Hlhv));
y.engine.torque = M_ice;
y.engine.power = P_ice;
y.engine.powerLimit = [(param.cPice(1)*w_ice^2 + param.cPice(2)*w_ice + param.cPice(3)); ....
                       (param.cPice(4)*w_ice^2 + param.cPice(5)*w_ice + param.cPice(6))];

y.exhaust.pressure = p_em;
y.exhaust.temperature = T_em;

y.turbine.speed = w_tc;
y.turbine.pressureRatio = Pi_t;
y.turbine.massflow = dot_m_t;
y.turbine.BSR = BSR;
y.turbine.BSRMax = param.BSR_max;
y.turbine.BSRMin = param.BSR_min;
y.turbine.efficiency = eta_tm;
y.turbine.power = P_t_eta_tm;

y.wastegate.control = u_wg;
y.wastegate.massflow = dot_m_wg;

y.turbocharger.speed = w_tc;

y.generator.power = P_gen;
y.generator.energy = E_gen;


end

function param = gensetParameters()
param.AFs= 14.7500;
param.A_t_eff= 7.8800e-04;
param.A_wg_eff= 3.1400e-04;
param.BSR_max= 1.1000;
param.BSR_min= 0.1500;
param.BSR_opt= 0.6220;
param.Hlhv= 42900000;
param.J_genset= 2.5000;
param.J_tc= 1.8695e-04;
param.P_ice_max= 205000;
param.Phi_opt= 0.0598;
param.Psi_max= 1.6380;
param.Q= [79.2671 0.7985 -1.4358];
param.R_a= 287;
param.R_c= 0.0328;
param.R_e= 286;
param.R_t= 0.0325;
param.T_amb= 298.4636;
param.T_amb_r= 298.8639;
param.T_c_b= 295.5297;
param.T_im= 292.0892;
param.V_D= 0.0067;
param.V_em= 0.0015;
param.V_is= 0.0143;
param.V_tot= 1;
param.cPice= [3.7461 623.4000 -21721 -3.2003 1.8788e+03 -62537];
param.c_eta_vol= [1.4393e-04 -0.0239 1.1339];
param.c_fr= [3.7591e-05 -0.0052 0.7196];
param.c_m= 2.0245;
param.c_mc_surge= [12.3530 0.6861];
param.c_t= [0.6859 2.3633];
param.c_wg= [0.6666 5.3517];
param.cp_a= 1011;
param.cp_e= 1332;
param.cv_a= 724;
param.dot_m_c_corr_max= [0.0419 0.3416 0.1918];
param.eta_c_max= 0.7915;
param.eta_ig_isl= [0.6650 0.7090 0.4450];
param.eta_sc= 1.2563;
param.eta_tm_max= 0.7075;
param.gamma_a= 1.3964;
param.gamma_cyl= 1.3500;
param.gamma_e= 1.2734;
param.gen2= [2.1584e-06 -6.2603e-04 1.0956 -1.5119e-06 8.0325e-04 0.8229 360.3522];
param.h_tot= 44.6091;
param.lambda_min= 1.2;
param.n_cyl= 6;
param.p_amb= 1.0111e+05;
param.p_c_b= 9.9510e+04;
param.p_es= 1.0587e+05;
param.r_c= 17.2000;
param.u_max= [150; 1; 300000];
param.u_min= [0; 0; 0];
param.w_tc_corr_opt= 0.5106;
end