classdef ast_cat < yop.ast_expression
    properties
        d
        args
    end
    methods
        function obj = ast_cat(d, varargin)
            obj@yop.ast_expression();
            obj.d = d;
            obj.args = varargin;
            tmp = cell(size(varargin));
            for k=1:length(varargin)
                tmp{k} = ones(size(varargin{k}));
            end
            obj.dim = size(cat(d, tmp{:}));
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            tmp = {isa_numeric(obj.args{1})};
            for k=2:length(obj.args)
                tmp = {tmp{:}, isa_numeric(obj.args{k})};
            end
            boolv = cat(obj.d, tmp{:});
        end
        
        function obj = set_pred(obj)
            for k=1:length(obj.args)
                add_pred(obj.args{k}, obj);
            end
        end
        
        function boolv = isa_der(obj)
            tmp = {isa_der(obj.args{1})};
            for k=2:length(obj.args)
                tmp = {tmp{:}, isa_der(obj.args{k})};
            end
            boolv = cat(obj.d, tmp{:});
        end
        
        function [bool, id] = isa_variable(obj)
            [bool, id] = isa_variable(obj.args{1});
            tmp_bool = {bool};
            tmp_id = {id};
            for k=2:length(obj.args)
                [bk, ik] = isa_variable(obj.args{k});
                tmp_bool = {tmp_bool{:}, bk};
                tmp_id = {tmp_id{:}, ik};
            end      
            bool = cat(obj.d, tmp_bool{:});
            id = cat(obj.d, tmp_id{:});
        end
        
        function [bool, tp] = isa_timepoint(obj)
            [bool, tp] = isa_timepoint(obj.args{1});
            tmp_bool = {bool};
            tmp_tp = {tp};
            for k=2:length(obj.args)
                [bk, tk] = isa_timepoint(obj.args{k});
                tmp_bool = {tmp_bool{:}, bk};
                tmp_tp = {tmp_tp{:}, tk};
            end
            bool = cat(obj.d, tmp_bool{:});
            tp = cat(obj.d, tmp_tp{:});
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.args{k});
            end
            value = cat(evaluate(obj.d), tmp{:});
        end
        
        function v = forward(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = value(obj.args{k});
            end
            obj.m_value = cat(value(obj.d), tmp{:});
            v = obj.m_value;
        end
        
        function draw(obj)
            % every vararg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['cat(dim, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            draw(obj.d);
            end_child(obj);
            for k=1:(length(obj.args)-1)
                begin_child(obj);
                draw(obj.args{k});
                end_child(obj);
            end
            last_child(obj);
            draw(obj.args{end});
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
            [topsort, n_elem, visited] = ...
                topological_sort(obj.d, visited, topsort, n_elem);
            
            for k=1:length(obj.args)
                [topsort, n_elem, visited] = topological_sort(...
                    obj.args{k}, ...
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