classdef settings < handle
    
    properties
        m_errors = true
        m_warnings = true
        m_cx_type = yop.settings.MX
    end
    
    properties (Constant, Hidden)
        MX = 1
        SX = 2
    end
    
    methods
        function obj = settings()
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
        
        function bool = errors(value)
            obj = yop.settings();
            if nargin == 1
                obj.m_errors = value;
            end
            bool = obj.m_errors;
        end
        
        function type = cx_type(value)
            obj = yop.settings();
            if nargin == 1
                obj.m_cx_type = value;
            end
            type = obj.m_cx_type;
        end
        
    end
end