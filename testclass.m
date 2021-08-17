classdef testclass < handle
    properties
    end
    methods
        function obj = testclass()
        end
        
        function varargout = subsref(obj, s)
            
            if length(s) > 1 || s(1).type ~= "()"
                [varargout{1:nargout}] = builtin('subsref',obj, s);
                return;
            end
            % From this point length(s) == 1 && s.type == "()"

            % Determine if subs contains only numeric elements or colon.
            % subs: {[s1], [s2], ..., [sN]}
            numeric = true;
            for k=1:length(s.subs)
                numeric = (isnumeric(s.subs{k}) || ischar(s.subs{k}) ) ...
                    && numeric;
            end
            
            if length(s.subs)==1 && isa(s.subs{1}, 'yop.ast_eq')
                % timed expression: obj(t==4)
                varargout{1} = yop.ast_timepoint(s.subs{1}, obj);
                
            else
                % Use built-in for any other expression
                [varargout{1:nargout}] = builtin('subsref',obj, s);
                
            end

        end
    end
end