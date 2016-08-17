classdef xsClass
    %XSCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        scatOrder
        total
        transport
        absorption
        nufission
        fission
        chi
        scatter
    end
    
    methods
        function obj = xsClass2( name, order )
            obj.name = name;
            obj.scatOrder = order;
        end
        
        function obj = xsClass( ngroups )
            obj.transport(1:ngroups,1) = 0.0;
        end
    end
    
end

