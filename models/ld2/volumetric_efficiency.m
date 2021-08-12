function eta_vol = volumetric_efficiency(p_im, N_ice, param)

eta_vol = param.c_vol(1)*sqrt(p_im) + ...
          param.c_vol(2)*sqrt(N_ice) + ...
          param.c_vol(3);
end

