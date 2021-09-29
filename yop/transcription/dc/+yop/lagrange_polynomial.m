classdef lagrange_polynomial < handle & matlab.mixin.Copyable
    % LAGRANGE_POLYNOMIAL A class for creating and using Lagrange
    % polynomials.
    %    The class implements basic functionality for interpolating a set
    %    of values using Lagrange polynomials. The Lagrange polynomials are
    %    calculated as:
    %
    %       L(x) = sum( l_j(x)*y_j )_{j=0}^{d}
    %       l_j(x) = prod( (x - x_r)/(x_j - x_r) )_{r=0, r/=j}^{d}.
    %
    %    where L is the Lagrange polynomial, l_j the basis polynimals,
    %    x the independent variable, x_i the sampling timepoints, y_j the
    %    sampled values, and d the polynomial degree.
    %
    % -- Properties --
    %    x : Row vector describing the sampling timepoints
    %
    %    y : Matrix describing the sampling values. Must have equally many
    %        columns as in 'x'. E.g. if the sampled signal is scalar-
    %        valued, 'y' is a row vector
    %
    %    l : Matrix containing the Lagrange basis polynomials.
    %
    % -- Methods --
    %    obj = lagrange_polynomial() : Constructor.
    %
    %    obj = calculate_basis(obj) : Calculates the basis polynomials.
    %
    %    values = evaluate(obj, t) : Evaluate the polynomial at time t.
    %
    %    polynomial=integrate(obj,constant_term):Integrates the polynomial.
    %
    %    polynomial = differentiate(obj) : Differentiate the polynomial.
    %
    %    deg = degree(obj) : Get the polynomial degree.
    %
    % -- Examples --
    %    % COPY INTO SCRIPT
    %    % Approximate the t^2 and t^3  at the specified timepoints using a
    %    % second order Lagrange polynomial.
    %    timepoints = [0,1,2,3]; % Sample points results in 3nd ord. polyn.
    %    analytical_values = @(t) [t.^2; t.^3];
    %    analytical_derivative = @(t) [2*t; 3*t.^2];
    %    analytical_integral = @(t) [t.^3/3; t.^4/4];
    %
    %    lp = yop.lagrange_polynomial(timepoints, analytical_values(timepoints));
    %
    %    t = 0:0.05:4;
    %    figure(1); hold on
    %    plot(t, analytical_values(t))
    %    plot(t, lp.evaluate(t), 'x')
    %    legend('t^2', 't^3', 'lp_1', 'lp_2')
    %    title('Polynomial approximation')
    %
    %    figure(2); hold on
    %    plot(t, analytical_derivative(t))
    %    plot(t, lp.differentiate.evaluate(t), 'x')
    %    legend('2*t', '3*t.^2', 'lp_1', 'lp_2')
    %    title('differentiation')
    %
    %    figure(3); hold on
    %    plot(t, analytical_integral(t))
    %    plot(t, lp.integrate.evaluate(t), 'x')
    %    legend('1/3 t.^3', '1/4 t.^4', 'lp_1', 'lp_2')
    %    title('integration')
    %
    % -- Details --
    %    For details, see:
    %    https://en.wikipedia.org/wiki/Lagrange_polynomial
    
    properties
        x  % Sample timepoints
        y  % Sample values
        l  % Polynomial basis
    end
    methods
        
        function obj = lagrange_polynomial(x, y)
            % LAGRANGE_POLYNOMIAL Constructor
            %    Constructs the polynomial from the sampling points x and
            %    the sampled values y.
            %
            % -- Syntax --
            %    obj = yop.lagrange_polynomial(x, y);
            %
            % -- Parameters --
            %    x : Sampled timepoints. Specified as a row vector.
            %
            %    y : Sampled values. Specified as a Matrix.
            %        Must have as many columns as 'x' do.
            %        Dimension of values are specified in the
            %        column direction.
            %
            % -- Examples --
            %    % Scalar values
            %    obj = yop.lagrange_polynomial([1 2 3], [1 4 9])
            %
            %    % Vector valued values
            %    obj = yop.lagrange_polynomial([1 2 3], [1 4 9; 1 8 27])
            
            assert(size(x,1)==1, ...
                '[Yop] Error: Argument x is not a row vector.');
            
            assert(size(x,2)==size(y,2), ['[Yop] Error: Argument y ', ...
                'does not have equally many columns as x.']);
            
            obj.x = x;
            obj.y = y;
            calculate_basis(obj);
        end
        
        function obj = calculate_basis(obj)
            % CALCULATE_BASIS Calculates the Lagrange polynomial basis
            %    For the given Lagrange polynomial:
            %
            %       L(x) = sum( l_j(x)*y_j )_{j=0}^{d}
            %
            %    calculates the basis l_j for all j according to:
            %
            %       l_j(x) = prod( (x - x_r)/(x_j - x_r) )_{r=0, r/=j}^{d}.
            %
            % -- Syntax --
            %    obj.calulate_basis()
            %    calculate_basis(obj)
            %
            % -- Parameters --
            %    obj : Handle to the Lagrange polynomial.
            
            obj.l = zeros(obj.degree+1, obj.degree+1);
            for j=1:obj.degree+1
                l_j = 1;
                for r=1:obj.degree+1
                    if j~=r
                        Pi_r = [1 -obj.x(r)] / (obj.x(j)-obj.x(r));
                        l_j = conv(l_j, Pi_r);
                    end
                end
                obj.l(j,:) = yop.lagrange_polynomial.filter0(l_j);
            end
        end
        
        function values = evaluate(obj, t)
            % EVALUATE Evaluates the polynomial at time t
            %    Evaluate the polynomial at time t, by evaluating the
            %    Lagrange polynomial:
            %
            %       L(t) = sum( l_j(t)*y_j )_{j=0}^{d}
            %
            %    If t is a vector, it evaluates the polynomial at all
            %    timepoints.
            %
            % -- Syntax --
            %    obj.evaluate(t)
            %    evaluate(obj, t)
            %
            % -- Parameters --
            %    obj : Handle to the Lagrange polynomial instance.
            %
            %    t   : Vector with the timepoints the polynomial should be
            %           evaluated at.
            %
            % -- Examples --
            %    lp.evaluate(1);
            %    lp.evaluate(0:0.1:1);
            %    evaluate(lp, 2);
            %    evaluate(lp, 1:10);
            
            values = [];
            for n=1:length(t)
                v = 0;
                for j=1:obj.degree+1
                    ll = yop.lagrange_polynomial.filter0( ...
                        polyval(obj.l(j,:), t(n)));
                    v = v + ll .* obj.y(:,j);
                end
                values = [values, v];
            end
        end
        
        function polynomial = integrate(obj, constant_term)
            % INTEGRATE Integrates the Lagrange polynomial
            %    Integrates the Lagrange polynomial with an optional
            %    constant term. Returns a new lagrange polynomial that
            %    is the integration of the input.
            %
            % -- Syntax --
            %    polynomial = integrate(obj, constant_term)
            %    polynomial = obj.integrate(constant_term)
            %
            % -- Parameters --
            %    obj           : Handle to the Lagrange polynomial to be
            %                    integrated.
            %
            %    polynomial    : Integration of the input polynomial
            %                    described as a new Lagrange polynomial.
            %
            % -- Arguments (Optional) --
            %    constant_term : Constant of integration, specified as a
            %                    numeric scalar. Defaults to 0.
            %
            % -- Examples --
            %    lp_int = lp.integrate(0)
            %    lp_int = integrate(lp, 0)
            
            if nargin == 1
                constant_term = 0;
            end
            
            [r, c] = size(obj.l);
            new_basis = zeros(r, c+1);
            for k=1:r
                new_basis(k,:) = polyint(obj.l(k,:), constant_term);
            end
            polynomial = copy(obj);
            polynomial.l = yop.lagrange_polynomial.filter0(new_basis);
        end
        
        function polynomial = differentiate(obj)
            % DIFFERENTIATE Differentiate the Lagrange polynomial
            %    Differentiates the Lagrange polynomial obj and returns a
            %    new lagrange polynomial that is the differentiated
            %    polynomial of the input.
            %
            % -- Syntax --
            %    polynomial = obj.differentiate()
            %    polynomial = differentiate(obj)
            %
            % -- Parameters --
            %    obj        : Handle to the Lagrange polynomial to be
            %                 differentiated.
            %
            %    polynomial : Differentiation of the input polynomial
            %                 described as a new Lagrange Polynomial.
            %
            % -- Examples --
            %    lp_der = lp.differentiate()
            %    lp_der = differentiate(lp)
            
            [r, c] = size(obj.l);
            new_basis = zeros(r, c-1);
            for k=1:r
                new_basis(k, :) = polyder(obj.l(k, :));
            end
            polynomial = copy(obj);
            polynomial.l = yop.lagrange_polynomial.filter0(new_basis);
        end
        
        function deg = degree(obj)
            % DEGREE Get the Lagrange polynomial degree
            %
            % -- Syntax --
            %     deg = obj.degree
            %     deg = degree(obj)
            %
            % -- Parameters --
            %    obj : Handle to the Lagrange polynomial
            %    deg : The degree of the polynomial
            %
            % -- Examples --
            %    d = lp.degree
            %    d = degree(lp)

            deg = size(obj.x, 2)-1;
        end
        
    end
    
    methods (Static)
        function A = filter0(A)
            % FILTER0 Remove numerical noice from values close to zero
            % 
            % -- Syntax --
            %     A = yop.lagrange_polynomial.filter0(A)
            %
            % -- Parameters --
            %    A : Values to be filtered
            %
            for n=1:numel(A)
                if yop.EQ(A(n), 0, 1e-12) % floating point comparison
                    A(n) = 0;
                end
            end
        end
    end
    
end