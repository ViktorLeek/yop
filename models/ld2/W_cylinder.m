function [ W_cyl ] = W_cylinder(p_im, T_im, N_ice, param )

W_cyl = param.eta_vol * param.V_D .* p_im .* N_ice ./ ...
        (param.n_r * param.R_a .* T_im * 60);


end

