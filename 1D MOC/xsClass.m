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
        scatter % Scatter from column into row
    end
    
    methods
        function obj = xsClass2( name, order )
            obj.name = name;
            obj.scatOrder = order;
        end
        
        function obj = xsClass1( ngroups )
            obj.transport(1:ngroups,1) = 0.0;
        end
        
        function obj = calcTXS( obj, transOpt )
            obj.total = obj.absorption + obj.fission + sum(obj.scatter,1);
            switch(transOpt)
                case('P0')
                    obj.transport = obj.total;
            end
        end
    end
    
end

