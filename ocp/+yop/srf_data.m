classdef srf_data < handle
    properties
        lt = {}
        gt = {}
        le = {}
        ge = {}
        eq = {}
        ne = {}
    end
    methods
        function add_lt(obj, r)
            obj.lt = {obj.lt{:}, r};
        end
        
        function add_gt(obj, r)
            obj.gt = {obj.gt{:}, r};
        end
        
        function add_le(obj, r)
            obj.le = {obj.le{:}, r};
        end
        
        function add_ge(obj, r)
            obj.ge = {obj.ge{:}, r};
        end
        
        function add_eq(obj, r)
            obj.eq = {obj.eq{:}, r};
        end
        
        function add_ne(obj, r)
            obj.ne = {obj.ne{:}, r};
        end
        
        function r = get_relations(obj)
            r = {obj.lt{:}, obj.gt{:}, obj.le{:}, obj.ge{:}, ...
                obj.eq{:}, obj.ne{:}};
        end
    end
end