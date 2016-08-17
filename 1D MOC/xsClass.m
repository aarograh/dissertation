classdef xsClass
    %XSCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        transport
    end
    
    methods
        function obj = xsClass( ngroups )
            obj.transport(1:ngroups,1) = 0.0;
        end
    end
    
end

