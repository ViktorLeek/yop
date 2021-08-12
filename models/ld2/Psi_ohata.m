function [Psi] = Psi_ohata(Pi, gamma)

Pi_choke = 1/(gamma+1);
Psi_choke = sqrt( (gamma+1)/(2*gamma) .* (1-Pi_choke) .* (Pi_choke + (gamma-1)./(gamma+1)) );
Psi_nonsat = sqrt( (gamma+1)/(2*gamma) .* (1-Pi) .* (Pi + (gamma-1)./(gamma+1)) );

sp = logistic_function(Pi, Pi_choke, 1, 80);
Psi = Psi_choke + sp.*(Psi_nonsat - Psi_choke);

% % 
% Pi(Pi < 1/(gamma+1)) = 1/(gamma+1);
% Psi = sqrt( (gamma+1)/(2*gamma) .* (1-Pi) .* (Pi + (gamma-1)/(gamma+1)) );

end
