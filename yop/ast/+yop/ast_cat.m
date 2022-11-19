classdef ast_cat < yop.ast_expression
    properties
        d
        args
        m_dval
    end
    methods
        function obj = ast_cat(d, varargin)
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
                cat(d, val{:}), ... value
                cat(d, num{:}), ... numval
                cat(d, tt0{:}), ... t0
                cat(d, ttf{:}), ... tf
                cat(d, der{:}), ... der
                cat(d, red{:}), ... reducible
                cat(d, typ{:}), ... type
                cat(d, tid{:}) ... typeid
                );
            obj.m_d = d;
            obj.m_args = varargin;
            if all(obj.m_type == yop.var_type.state)
                dval = c0;
                for k=1:length(varargin)
                    dval{k} = varargin{k}.m_dval;
                end
                obj.m_dval = cat(d, dval{:});
            end
        end
        
        function d = der(obj)
            if isempty(obj.m_dval)
                d = der@yop.ast_expression(obj);
            else
                d = obj.m_dval;
            end
        end
        
        function ast(obj)
            % every vararg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.m_args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['cat(dim, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.m_d);
            end_child(obj);
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
            if ~isempty( find(visited == obj.m_id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.m_id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_d, visited, topsort, n_elem);
            
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