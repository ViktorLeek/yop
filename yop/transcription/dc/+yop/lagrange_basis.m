classdef lagrange_basis < handle & matlab.mixin.Copyable
    properties
        tau % including 0!
        l
    end
    methods
        function obj = lagrange_basis(tau)
            obj.tau = tau;
            
            % Compute basis
            obj.l = zeros(obj.n_points, obj.n_points);
            for j=1:obj.n_points
                l_j = 1;
                for r=1:obj.n_points
                    if j~=r
                        Pi_r = [1 -obj.tau(r)] / (obj.tau(j) - obj.tau(r));
                        l_j = conv(l_j, Pi_r);
                    end
                end
                obj.l(j,:) = l_j;
            end
        end
        
        function D = eval_basis(obj, t)
            D = zeros(obj.degree+1, length(t));
            for n=1:length(t)
                for j=1:obj.degree+1
                    D(j,n) = polyval(obj.l(j,:), t(n));
                end
            end
            D = langrange_polynomial.filter_zeros(D);
        end
        
        function polynomial = integrate(obj, constant_term)
            if nargin == 1
                constant_term = 0;
            end
            
            [r, c] = size(obj.l);
            new_basis = zeros(r, c+1);
            for k=1:r
                new_basis(k,:) = polyint(obj.l(k,:), constant_term);
            end
            polynomial = copy(obj);
            polynomial.l = langrange_polynomial.filter_zeros(new_basis);
        end
        
        function polynomial = differentiate(obj)
            [r, c] = size(obj.l);
            new_basis = zeros(r, c-1);
            for k=1:r
                new_basis(k, :) = polyder(obj.l(k, :));
            end
            polynomial = copy(obj);
            polynomial.l = langrange_polynomial.filter_zeros(new_basis);
        end
        
        
        function d = n_points(obj)
            d = length(obj.tau);
        end
    end
    
    methods (Static)
        function A = filter_zeros(A)
            % Filter out numerical noice and set values sufficiently close
            % to zero to zero.
            for n=1:numel(A)
                if yop.EQ(A(n), 0, 1e-12)
                    A(n) = 0;
                end
            end
        end
        
        function value = evaluate(basis, params)
            
        end
    end
end