function bool = EQ(a, b, tol)
if nargin == 2
    tol = 1e-9;
end
bool = abs(a-b) <= tol;
end