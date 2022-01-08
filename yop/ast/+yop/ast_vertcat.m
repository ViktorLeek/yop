classdef ast_vertcat < yop.ast_expression
    properties
        m_args
    end
    methods
        function obj = ast_vertcat(varargin)
            c0 = cell(size(varargin));
            val=c0; num=c0; tt0=c0; ttf=c0; der=c0; red=c0; typ=c0; tid=c0;
            for k=1:length(varargin)
                vk = varargin{k};
                val{k} = value(vk);
                num{k} = numval(vk);
                tt0{k} = get_t0(vk);
                ttf{k} = get_tf(vk);
                der{k} = get_der(vk);
                red{k} = isa_reducible(vk);
                [typk, tidk] = Type(vk);
                typ{k} = typk;
                tid{k} = tidk;
            end
            obj@yop.ast_expression( ...
                vertcat(val{:}), ... value
                vertcat(num{:}), ... numval
                vertcat(tt0{:}), ... t0
                vertcat(ttf{:}), ... tf
                vertcat(der{:}), ... der
                vertcat(red{:}), ... reducible
                vertcat(typ{:}), ... type
                vertcat(tid{:}) ... typeid
                );
            obj.m_args = varargin;
        end
        
        function ast(obj)
            % every arg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.m_args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['vertcat(', str(1:end-2), ')\n']);
            
            for k=1:(length(obj.m_args)-1)
                begin_child(obj);
                ast(obj.m_args{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.m_args{end});
            end_child(obj);
            
        end
        
        function [topsort, n_elem, visited] = ...
                topological_sort(obj, visited, topsort, n_elem)
            % Topological sort of expression graph by a dfs.
            
            switch nargin
                case 1
                    % Start new sort: topological_sort(obj)
                    visited = [];
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                case 2
                    % Semi-warm start: topological_sort(obj, visited)
                    % In semi-warm start 'visited' is already provided, but 
                    % no elements are sorted. This is for instance useful 
                    % for finding all variables in a number of expressions 
                    % that are suspected to contain common subexpressions.
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                otherwise
                    % Pass
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            for k=1:length(obj.m_args)
                [topsort, n_elem, visited] = topological_sort(...
                    obj.m_args{k}, ...
                    visited, ...
                    topsort, ...
                    n_elem ...
                    );
            end
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end