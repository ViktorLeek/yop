function y = logistic_function(x, x0, L, k)
y = L./(1+exp( -k*(x-x0) ));
end