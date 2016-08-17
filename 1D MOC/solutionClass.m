classdef solutionClass
    %SOLUTIONCLASS Stores angular and scalar flux solution variables
    %   This is mainly to make it easier to pass information around.
    %   It will also be used to check convergence for eigenvalue
    %   problems.
    
    properties
        BCond
        keff % First index is current, second is previous iteration
        angflux
        angflux_cellavg
        scalflux % First index is current, second is previous iteration
    end
    
    methods
        function obj = solutionClass( ncells,npol,ngroups,BCond )
            %SOLUTIONCLASS Constructor for solutionClass
            %   ncells  - The number of cells in the problem
            %   npol    - The number of polar angles in the problem
            %   ngroups - The number of energy groups in the problem
            %   BCond   - The boundary condition for the angular flux
            
            obj.keff(1:2) = 1.0;
            obj.angflux(1:ncells+1,1:npol,1:2,1:ngroups) = 0.0;
            obj.scalflux(1:ncells,1:ngroups,1:2) = 0.0;
            obj.BCond = BCond;
        end
        
        function obj = update( obj )
            %UPDATEBC Updates the solution to prepare for the next iteration
            %   Resets the rest of the angular flux array to 0.0
            %   Applies the appropriate boundary conditions to the angular flux
            %   Saves the scalar flux before zeroing it
            %   TODO: update eigenvalue
            
            obj.angflux(2:end,:,1,:) = 0.0;
            obj.angflux(1:end-1,:,2,:) = 0.0;
            obj.angflux_cellavg = 0.0;
            obj.scalflux(:,:,2) = obj.scalflux(:,:,1);
            obj.keff(2) = obj.keff(1);
            
            % Set angular flux BC
            if ischar(obj.BCond(1))
                if strcmp(obj.BCond(1),'vacuum')
                    obj.angflux(1,:,1,:) = 0.0;
                end
            end
            if ischar(obj.BCond(2))
                if strcmp(obj.BCond(2),'vacuum')
                    obj.angflux(end,:,2,:) = 0.0;
                end
            end
            
            %TODO: update eigenvalue
        end
    end
    
end

