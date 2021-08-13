function [T_em, T_e] = T_em_model(p_im, p_em, T_im, W_f, W_cyl, param)

% Cylinder out temperature
T_e = T_cyl_out( p_im, p_em, T_im, W_f, W_cyl, param );

% Exhaust manifold temperature
T_em = param.T_amb + (T_e - param.T_amb).*exp( -(param.c_emh) ./ (param.c_pe.*(W_f + W_cyl)) ) ;

end