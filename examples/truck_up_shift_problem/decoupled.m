function [dx, y, c] = decoupled(x, u)

upshift_param;
w_tr=x(5); theta=x(6); w_w=x(7);

% MVEM dynamics and constratins, change intertia depending on phase
[dx, y, c] = mvem(x, u, 0, J_genset_phase(2));

% Driveline dynamcis
T_ds = k_ds * theta + c_ds * (w_tr/i_f - w_w);
F_a  = 0.5 * rho_air * c_d * A_f * w_w^2 * r_w^2;

dw_tr  = (-c_tr*w_tr - T_ds/i_f)/J_tr1;
dtheta = (w_tr/i_f - w_w);
dw_w   = (T_ds - F_a*r_w - F_r*r_w)/J_veh;
dx     = [dx; dw_tr; dtheta; dw_w; u(3)];

% Aux calcs
y.dtheta = dtheta;
y.dw_tr  = dw_tr;
y.M_ds   = T_ds;

end