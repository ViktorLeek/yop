classdef settings < handle
    properties
        m_warnings = true
    end
    methods
        function obj = settings(obj)
            persistent OBJ
            if isempty(OBJ)
                OBJ = obj;
            else
                obj = OBJ;
            end
        end
        
    end
    
    methods (Static)
        function bool = warnings(value)
            obj = yop.settings();
            if nargin == 1
                obj.m_warnings = value;
            end
            bool = obj.m_warnings;
        end
    end
end