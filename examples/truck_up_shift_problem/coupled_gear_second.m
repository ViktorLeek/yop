function [dx, y, c] = coupled_gear_second(x, u)

upshift_param;
w_e=x(1); theta=x(6); w_w=x(7);

% Air and drag resistnace
F_a = 0.5 * rho_air * c_d * A_f * w_w^2 * r_w^2;

% Gear ratio
i_tot = i_t(2) * i_f;

% Drive shaft torque
M_ds = k_ds*theta + c_ds * (w_e/i_tot - w_w);

% MVEM dynamics and constratins, change intertia depending on phase
[dx, y, c] = mvem(x, u, M_ds/i_tot, J_genset_phase(3));

% Dynamics not part of MVEM
dw_tr  = dx(1)/i_t(2);
dtheta = (w_e/i_tot - w_w);
dw_w   = (M_ds - F_a*r_w - F_r*r_w)/J_veh;

% Dynamics
dx = [dx; dw_tr; dtheta; dw_w; u(3)];

% Signals
y.dtheta = dtheta;
y.dw_tr = dw_tr;
y.M_ds = M_ds;

end