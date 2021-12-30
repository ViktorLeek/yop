clear

yopvar states: x size: [3,3] scaling: [1, 2, 3] offset: [-1, -2, -3]
yopvar controls: u1 u2 degree: [0, 2]

vars = [yop.ocp_state_control(u1), yop.ocp_state_control(u2)];
aug_states  = yop.ocp_state_control.empty(1,0);
aug_eqs = {};