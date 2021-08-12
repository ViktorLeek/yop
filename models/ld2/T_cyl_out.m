function [ T_cyl_out ] = T_cyl_out( p_im, p_em, T_im, W_f, W_a, param )
% [ T_cyl_out ] = T_cyl_out( p_im, p_em, T_im, W_f, W_a, param )

q_in = W_f .* param.q_HV ./ (W_f + W_a);
T_cyl_out = param.eta_sc .* (p_em./p_im).^(1 - 1/param.gamma_air) .* param.r_c^(1 - param.gamma_air) .*...
    (q_in/param.c_pa + T_im .* param.r_c^(param.gamma_air - 1));


end

