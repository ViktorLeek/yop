classdef ast_subsasgn < yop.ast_expression
    
    properties
        node
        s
        b
    end
    
    methods
        function obj = ast_subsasgn(node, s, b)
            
            % Casadi is a bit picky, so logical indices must be converted
            % to indices.
            for k=1:length(s.subs)
                if islogical(s.subs{k})
                    s.subs{k} = find(s.subs{k});
                end
            end
            obj@yop.ast_expression();
            obj.node = node;
            obj.s = s;
            obj.b = b;
            obj.dim = size(node);
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used, or an
            % instance variable to remember the answer (assumes graph is
            % static)
            boolv = isa_numeric(obj.node);
            idx = get_indices(obj);
            boolv(idx) = isa_numeric(obj.b);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = is_transcription_invariant(obj.node);
            idx = get_indices(obj);
            boolv(idx) = is_transcription_invariant(obj.b);
        end
        
        function obj = set_pred(obj)
            add_pred(obj.node, obj);
            add_pred(obj.b, obj);
        end
        
        function [bool, id] = isa_der(obj)
            [bool, id] = isa_der(obj.node);
            [boolb, idb] = isa_der(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_variable(obj)
            [bool, id] = isa_variable(obj.node);
            [boolb, idb] = isa_variable(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_state(obj)
            [bool, id] = isa_state(obj.node);
            [boolb, idb] = isa_state(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_independent(obj)
            [bool, id] = isa_independent(obj.node);
            [boolb, idb] = isa_independent(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_parameter(obj)
            [bool, id] = isa_parameter(obj.node);
            [boolb, idb] = isa_parameter(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_control(obj)
            [bool, id] = isa_control(obj.node);
            [boolb, idb] = isa_control(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, id] = isa_algebraic(obj)
            [bool, id] = isa_algebraic(obj.node);
            [boolb, idb] = isa_algebraic(obj.b);
            idx = get_indices(obj);
            bool(idx) = boolb;
            id(idx) = idb;
        end
        
        function [bool, tp] = isa_timepoint(obj)
            [bool, tp] = isa_timepoint(obj.node);
            idx = get_indices(obj);
            [boolb, tpb] = isa_timepoint(obj.b);
            bool(idx) = boolb;
            tp(idx) = tpb;
        end
        
        function value = evaluate(obj)
            % Subsref are only created if indices 's' are numerics, and
            % so they can be passed as they are.
            value = subsasgn(evaluate(obj.node), obj.s, evaluate(obj.b));
        end
        
        function v = forward(obj)
            obj.m_value = subsasgn(value(obj.node), obj.s, value(obj.b));
            v = obj.m_value;
        end
        
        function idx = get_indices(obj)
            % Return the indices this subsasgn node refers to
            node_val = obj.node.m_value;
            sz = size(obj.node);
            obj.node.m_value = reshape(1:prod(sz), sz);
            idx = subsref(value(obj.node), obj.s); % eval as subsasgn
            idx = idx(:);
            obj.node.m_value = node_val;
        end
        
        function draw(obj)
            fprintf('subsasgn(node, s, b)\n');
            
            begin_child(obj);
            draw(obj.node);
            end_child(obj);
            
            % Subs are numeric
            str = [];
            for k=1:length(obj.s.subs)
                str = [str, '[', num2str(obj.s.subs{k}), '], '];
            end
            begin_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
            end_child(obj);
            
            last_child(obj);
            draw(obj.b);
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
                topological_sort(obj.node, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.s, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.b, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end