function uf = smoke_limiter(u, u_fm, u_min, u_max)
% uf - Fuel injection
% u - Controller output
% u_fm - fuel max
% u_min - Actuator minimum value
% u_max - Actuator maximum value
uf = min(u, u_fm);  % Dynamic limit based on air flow
uf = max(uf, u_min);  % Static limit
uf = min(uf, u_max);  % Static limit
end