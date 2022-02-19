classdef ast_relation < yop.ast_binary_expression
    properties
        m_hard = false
    end
    
    methods
        function obj = ast_relation(value, numval, lhs, rhs)
            sz = size(value);
            t0_l = get_t0(lhs);
            t0_r = get_t0(rhs);
            tf_l = get_tf(lhs);
            tf_r = get_tf(rhs);
            t0 = max([t0_l(:); t0_r(:)]) * ones(sz);
            tf = min([tf_l(:); tf_r(:)]) * ones(sz);
            isder = false;
            isreducible = isa_reducible(lhs) & isa_reducible(rhs);
            type = zeros(size(sz));
            typeid = zeros(size(sz));
            obj@yop.ast_binary_expression(value, numval, t0, tf, isder, ...
                isreducible, type, typeid, lhs, rhs);
        end
        
        function obj = hard(obj)
            obj.m_hard = true;
        end
        
        function bool = is_hard(obj)
            bool = obj.m_hard;
        end
        
        function sz = size(obj, varargin)
            %             sz = [1, 1];
            sz = size(obj.m_value, varargin{:});
        end
        
        function rels = get_relations(obj)
            l = get_relations(obj.m_lhs);
            r = get_relations(obj.m_rhs);
            rels = {obj, l{:}, r{:}};
            rels = rels(~cellfun('isempty', rels));
        end
        
        function r = rmost(obj)
            r = rmost(obj.m_rhs);
        end
        
        function l = lmost(obj)
            l = lmost(obj.m_lhs);
        end
        
    end
end