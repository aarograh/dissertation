classdef xsClass
    %XSCLASS Class which hold cross-section data
    %   This class holds the multi-group XS data for
    %   a single type of material.
    
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
        function obj = xsClass( name, order )
            %XSCLASS Consructor for xsClass
            %   name  - name of the XS set
            %   order - Scattering order for the xsSet
            
            obj.name = name;
            obj.scatOrder = order;
            
        end
        
        function obj = calcTXS( obj, transOpt )
            %CALCTXS Calculates transport and total cross-sections
            %   obj - the XS set object
            %   transOpt - Transport correction option.  Currently only values
            %              of P0 are allowed.
            
            obj.total = obj.absorption + sum(obj.scatter,1);
            switch(transOpt)
                case('P0')
                    obj.transport = obj.total;
            end
            
        end
    end
    
end

