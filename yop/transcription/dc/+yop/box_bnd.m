function [w_lb, w_ub] = box_bnd(N, d, ocp)
t0_lb = ocp.t0_lb;
t0_ub = ocp.t0_ub;
tf_lb = ocp.tf_lb;
tf_ub = ocp.tf_ub;

reps = N*(d + 1) + 1;
x_lb = repmat(ocp.x_lb, reps, 1);
x_ub = repmat(ocp.x_ub, reps, 1);
x_lb(1 : ocp.n_x) = ocp.x0_lb;
x_ub(1 : ocp.n_x) = ocp.x0_ub;
x_lb(end - ocp.n_x + 1 : end) = ocp.xf_lb;
x_ub(end - ocp.n_x + 1 : end) = ocp.xf_ub;

u_lb = repmat(ocp.u_lb, N, 1);
u_ub = repmat(ocp.u_ub, N, 1);
u_lb(1 : ocp.n_u) = ocp.u0_lb;
u_ub(1 : ocp.n_u) = ocp.u0_ub;
u_lb(end - ocp.n_u + 1 : end) = ocp.uf_lb;
u_ub(end - ocp.n_u + 1 : end) = ocp.uf_ub;

p_lb = ocp.p_lb;
p_ub = ocp.p_ub;

w_ub = vertcat(t0_ub, tf_ub, x_ub, u_ub, p_ub);
w_lb = vertcat(t0_lb, tf_lb, x_lb, u_lb, p_lb);
end