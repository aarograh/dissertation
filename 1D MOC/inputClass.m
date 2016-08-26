classdef inputClass
    %INPUTCLASS Container class to hold input data
    %   Contains all data for pin geometry and materials,
    %   XS Library filename, and ray tracing and meshing
    %   parameters
    
    properties
        pitch
        pinmap
        diag
        pinmats
        radii
        pinmesh
        npol
        xsfilename
        scattype
        BCond
        nouters
        conv_crit=[0, 0]
        verbose
    end
    
    methods
        function obj = inputClass( )
            %INPUTCLASS Creates an empty inputClass
            %   No attributes are initialized.  This
            %   simply creates the class and leaves
            %   filling the attributes to the user.
            
        end
    end
    
end

