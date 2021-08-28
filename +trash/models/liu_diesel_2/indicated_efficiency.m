function [eta_ig, eta_otto, gamma_cyl, eta_cal, c_cal_1, c_cal_2] = indicated_efficiency(phi, N_ice, u_f, param)

% Calc gamma
gamma_cyl = param.c_gamma(1) + param.c_gamma(2)*phi + param.c_gamma(3)*phi.^2;

% Calc operating point dependent engine calibration factor
c_cal_1 = param.c_inj(2) + param.c_inj(3)*(N_ice/1000);
c_cal_2 = param.c_inj(4) + param.c_inj(5)*(N_ice/1000) + param.c_inj(6)*(N_ice/1000).^2;

eta_cal =  c_cal_2.*( u_f/100 - c_cal_1 ).^2 + param.c_inj(1);

% Indicated efficiency
eta_otto = 1 - 1./( param.r_c.^(gamma_cyl - 1));
eta_ig = eta_cal.*eta_otto;

end