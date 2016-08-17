classdef quadratureClass
    %QUADRATURECLASS Contains quadrature information
    %   Contains the angles, their cosines, and weights
    %   for a polar quadrature for 1D MOC
    
    properties
        npol
        angles
        cosines
        weights
    end
    
    methods
        function obj = quadratureClass(npol)
            %QUADRATURECLASS Sets up a quadrature class
            %   npol - Number of polar angles.  Currently, only values of
            %          1 are accepted.
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