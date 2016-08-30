classdef cmfdClass < handle
    %CMFDCLASS CMFD Accelerator class
    %   This class contains the data and methods required to perform
    %   1D CMFD acceleration
    
    properties
        dtils
        dhats
        A
        x
        b
        keff
        xsLib
        mesh
        ncells
        regPerCell
        firstIteration = true
    end
    
    methods
        function obj = cmfdClass( input, mesh, xsLib)
            %CMFDCLASS Constructor for cmfdClass object
            %   input - The inputClass object from which to initialize
            %   mesh  - The meshClass object to solve
            %   xsLib - The xsLibraryClass object with XS data
            
            npins = length(input.pinmap);
            nmats = size(input.pinmats,2);
            ngroups = xsLib.ngroups;
            obj.xsLib = xsLib;
            obj.mesh = mesh;
            obj.ncells = 0;
            for ipin=1:npins
                for imat=nmats:-1:1
                    if input.pinmats(input.pinmap(ipin),imat)
                        break
                    end
                end
                obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),imat:-1:1);
                obj.ncells = obj.ncells + imat;
                obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),1:imat);
                obj.ncells = obj.ncells + imat;
            end
            
            obj.dtils(1:obj.ncells+1,1:ngroups) = 0.0;
            obj.dhats(1:obj.ncells+1,1:ngroups) = 0.0;
            obj.A(1:obj.ncells*ngroups,1:obj.ncells*ngroups) = 0.0;
            obj.x(1:obj.ncells*ngroups,1) = 0.0;
            obj.b(1:obj.ncells*ngroups,1) = 0.0;
            
        end
        
        function obj = setup( obj, solution, mesh )
            %SETUP Sets up matrix and source for CMFD problem
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
        end
        
        function solution = solve( obj, solution, mesh )
            %SOLVE Performs power iteration to solve the CMFD system
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            obj.keff(1:2) = solution.keff(1);
            converged = 0;
            
            obj.homogenize(solution, mesh);
            while ~converged
                obj.step();
                converged = 1;
            end
            obj.project(solution, mesh);
            
            solution.keff(2) = solution.keff(1);
            solution.keff(1) = obj.keff(1);
            
        end
        
        function obj = step( obj )
            %STEP Performs a single CMFD iteration
            %   obj      - The cmfdClass object to set up
            
            
        end
        
        function obj = updateEig( obj )
            %UPDATEEIG Performs update of CMFD eigenvalue
            %   obj - cmfdClass object to update
            
        end
        
        function obj = project( obj, solution, mesh )
            %PROJECT Projects the CMFD solution onto the MOC mesh
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
        end
    end
    
end

