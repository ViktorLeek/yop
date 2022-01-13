function dw = calcWeight(v_cas, H, throttle, FFModel)

% coefficients for maximum throttle setting
am1 = -1.5374e-13;
am2 = +2.3875e-11;
am3 = -5.7989e-10;
bm1 = +2.1483e-07;
bm2 = -3.2731e-07;
bm3 = -9.1353e-07;
cm1 = -5.8602e-06;
cm2 = +9.7605e-04;
cm3 = +1.1389e-01; 

% coefficients for minimum throttle setting (i for idle)
ai1 = +4.2242e-15;        
ai2 = -9.0868e-13;       
ai3 = +1.6801e-10;    
bi1 = -4.8122e-11;   
bi2 = +1.0818e-08;     
bi3 = -2.7199e-06;    
ci1 = +2.5788e-08;     
ci2 = -6.8894e-06;    
ci3 = +2.5804e-02;        
  
% Maximum fuel flow
am = am1*v_cas^2 + am2*v_cas + am3;
bm = bm1*v_cas^2 + bm2*v_cas + bm3;
cm = cm1*v_cas^2 + cm2*v_cas + cm3;
Fmax = am*H^2 + bm*H + cm; % maximum fuel flow    [kg/s]

% Minimum fuel flow (idle throttle setting)
ai = ai1*v_cas^2 + ai2*v_cas + ai3;
bi = bi1*v_cas^2 + bi2*v_cas + bi3;
ci = ci1*v_cas^2 + ci2*v_cas + ci3;
F_idle = ai*H^2 + bi*H + ci; % minimum fuel flow [kg/s]

% Interpolation
FF = F_idle + (Fmax - F_idle)*ppval(FFModel, throttle); % Current fuel flow [kg/s]

% Calculation of weight
dw = -2*FF*9.81;

