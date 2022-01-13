classdef (InferiorClasses = {?yop.ocp}) multiphase < handle
    properties
        m_nlp
        m_ocps % in sequence
    end
    methods
        function obj = multiphase(varargin)
            K = length(varargin);
            for k=1:K-1
                
            end
            obj.m_nlp = varargin{K};
        end
    end
end