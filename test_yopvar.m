clear
% yopvar -time t -W 1 -os 0 -time0 t0 -w 2 -os -1 -timef tf -weight 3 -os -2
% yopvar -state x1 x2 x3 -w [1, 2, 3] -os [-1 -2 -3] -state x4 x5 -weight [4, 5] -os [-4, -5]
%yopvar -ctrl u -w 2 -os 3 -deg 4 -ctrl u1 u2 -w 2 -os 3 -deg [0,1] -state x1 x2 x3 -w [1, 2, 3] -os [-1 -2 -3] -state x4 x5 -weight 4 -os -5
% yopvar -state x1 x2 x3 -w [1, 2, 3] -os [-1 -2 -3] -state x4 x5 -weight 4 -os -5