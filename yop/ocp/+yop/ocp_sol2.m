classdef ocp_sol2 < handle
    properties
        t0
        tf
        t
        x
        z
        u
        p
    end
    methods
        function obj = ocp_sol2(t0, tf, t, x, z, u, p, N, tau)
                        
        end
        
        function [t,x,z,u,p] = structure_data(obj, w, d)
            
        end
        
        function ip = interpolating_poly(obj, t, v, N, tau)
            ip = yop.ocp_poly.empty(1,0);
            for n=1:N
                %                 ip(end+1) = yop.ocp_poly(
            end
        end
    end
end