classdef quadrature
    %QUADRATURE Summary of this class goes here
    %   Contains polar quadrature data for 1D MOC
    
    properties
        npol
        angles
        cosines
        weights
    end
    
    methods
        function obj = quadrature(npol)
            obj.npol = npol;
            if npol == 1
                obj.angles = pi/4;
                obj.weights = 2*pi;
                obj.cosines = cos(obj.angles);
            else
                fprintf('ERROR: Number of polar angles must be 1.')
                exit(npol)
            end
        end
    end
    
end