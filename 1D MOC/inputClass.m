classdef inputClass < handle
    %INPUTCLASS Container class to hold input data
    %   Contains all data for pin geometry and materials,
    %   XS Library filename, and ray tracing and meshing
    %   parameters
    
    properties
        % The pin pitch in the cylindrical model
        pitch
        % The map of pin IDs for the model
        pinmap
        % The list of materials for each pin
        %   The number of rows should be equal to the highest pin ID used in pinmap
        %   The number of columns should be equal to the largest number of radial regions in any pin
        pinmats
        % The list of radii for each pin
        %   The number of rows should be the same as pinmats
        %   The number of columns should be one less than pinmats
        radii
        % The number of radial subdivisions for each material region
        %   The shape should be the same as pinmats
        pinmesh
        % Number of polar angles to use in the Gaussian quadrature
        npol
        % Name of cross section library file
        xsfilename
        % Boundary condition.  Must be a vector of strings with length 2.  Acceptable options are 'vaccum' and
        %   'refelcting'.  The first entry is for the left boundary, while the second is for the right.
        BCond
        % Maximum number of outer iterations to solve
        nouters
        % Convergence criteria.  Must be vector of doubles with length 2
        %   First element is for change in k-eff
        %   Second element is for maximum change in fission source
        conv_crit=[0, 0]
        % Option to toggle the amount of output produced by the code
        verbose=false
        % The number of mixtures to be produced by the code.  These mixtures are volume homogenizations of
        %   materials specified in the cross section library file
        nmixtures=0
        % The material IDs to be mixed
        %   The number of rows is equal to nmixtures
        %   Each row is a list of material IDs to be mixed.  The number of columns is equal to the maximum
        %   number of materials being mixed.  All 0s ending a row will be truncated.
        mixtures
        % The volume fractions used to mix the materials
        %   The shape is the same as the shape of mixtures.  Each element of mixvols gives the volume fraction
        %   for the material in the same position in mixtures.
        mixvols
        % Toggles subray MOC on or off
        %   0: No subray MOC.  Just performs a single 2D MOC sweep using homogenized cross sections where
        %   necessary.
        %   1: Performs separate sweeps during each iteration, then combines the solutions afterwards using
        %   volume fractions.  This is the same as using subray MOC everywhere.
        %   2: Performs subray MOC in the regions determined by npinSubTrack.
        subray=0
        % Determines how many neighboring pins will use subray MOC if subray is set to 2.  A value of 0 will
        %   use subray MOC only in the pin cells that have homogenized cross sections.  Larger values add pin
        %   cells to each side to resolve more effects near the partially inserted rod.  Defaults to a very
        %   large value to use subray MOC everywhere (which corresponds to subray=1).
        npinSubTrack=10000000
        % The following 4 options are unused, but left for future work and legacy input support
        diag % Not supported
        cmfd % CMFD not supported
        voleql=false % Ignored
        scattype % Only isotropic scattering is currently supported
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

