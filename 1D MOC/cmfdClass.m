classdef cmfdClass < handle
    %CMFDCLASS CMFD Accelerator class
    %   This class contains the data and methods required to perform
    %   1D CMFD acceleration
    
    properties
        xstr
        xst
        xsrm
        xssc
        xsnf
        xsch
        dtils
        dhats
        A
        x
        b
        keff
        xsLib
        mesh
        ngroups
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
            
            obj.ncells = length(input.pinmap);
            obj.ngroups = xsLib.ngroups;
            obj.xsLib = xsLib;
            obj.mesh = mesh;
            for ipin=1:obj.ncells
                obj.regPerCell(ipin) = sum(input.pinmesh(input.pinmap(ipin),:),2);
            end
%             for ipin=1:npins
%                 for imat=nmats:-1:1
%                     if input.pinmats(input.pinmap(ipin),imat)
%                         break
%                     end
%                 end
%                 obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),imat:-1:1);
%                 obj.ncells = obj.ncells + imat;
%                 obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),1:imat);
%                 obj.ncells = obj.ncells + imat;
%             end
            
            obj.xstr(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xst(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xsrm(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xssc(1:obj.ncells,1:obj.ngroups,1:obj.ngroups) = 0.0;
            obj.xsnf(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xsch(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.dtils(1:obj.ncells+1,1:obj.ngroups) = 0.0;
            obj.dhats(1:obj.ncells+1,1:obj.ngroups) = 0.0;
            obj.A(1:obj.ncells*obj.ngroups,1:obj.ncells*obj.ngroups) = 0.0;
            obj.x(1:obj.ncells*obj.ngroups,1) = 0.0;
            obj.b(1:obj.ncells*obj.ngroups,1) = 0.0;
            
        end
        
        function solution = solve( obj, solution, mesh )
            %SOLVE Performs power iteration to solve the CMFD system
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            obj.keff(1:2) = solution.keff(1);
            converged = 0;
            
            obj.homogenize(solution, mesh);
            obj.setup();
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
        
        function obj = homogenize( obj, solution, mesh )
            %HOMOGENIZE Homogenizes the cmfd mesh
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            obj.xstr(:) = 0.0;
            obj.xst(:) = 0.0;
            obj.xsrm(:) = 0.0;
            obj.xssc(:) = 0.0;
            obj.xsnf(:) = 0.0;
            obj.xsch(:) = 0.0;
            
            for g=1:obj.ngroups
                icell = 0;
                for i=1:obj.ncells
                    volsum = 0.0;
                    flxvolsum = 0.0;
                    fisssrcsum = 0.0;
                    for j=1:obj.regPerCell(i)
                        icell = icell + 1;
                        imat = mesh.materials(icell);
                        dx = (mesh.fsredges(icell+1)-mesh.fsredges(icell));
                        volsum = volsum + dx;
                        flxvol = solution.scalflux(icell,g,1)*dx;
                        flxvolpsi = fisssrc(icell,1)*dx;
                        flxvolsum = flxvolsum + flxvol;
                        fisssrcsum = fisssrcsum + flxvolpsi;
                        obj.xstr(i,g) = obj.xstr(i,g) + obj.xsLib.xsSets(imat).transport(g)*flxvol;
                        obj.xst(i,g) = obj.xst(i,g) + obj.xsLib.xsSets(imat).transport(g)*flxvol;
                        for g2=1:obj.ngroups
                            obj.xssc(i,g2,g) = obj.xssc(i,g2,g) + ...
                                obj.xsLib.xsSets(imat).scatter(g2,g)*flxvol;
                        end
                        obj.xsnf(i,g) = obj.xsnf(i,g) + obj.xsLib.xsSets(imat).nufission(g)*flxvol;
                        obj.xsch(i,g) = obj.xsch(i,g) + ...
                            flxvolpsi*obj.xsLib.xsSets(imat).chi(g);
                    end
                    obj.xstr(i,g) = obj.xstr(i,g)/flxvolsum;
                    obj.xst(i,g) = obj.xst(i,g)/flxvolsum;
                    obj.xsnf(i,g) = obj.xsnf(i,g)/flxvolsum;
                    obj.xsch(i,g) = obj.xsch(i,g)/fisssrcsum;
                    for g2=1:obj.ngroups
                        obj.xssc(i,g2,g) = obj.xssc(i,g2,g)/flxvolsum;
                    end
                    obj.xsrm(i,g) = obj.xst(i,g) - obj.xssc(i,g,g);
                end
            end
            
        end
        
        function obj = setup( obj )
            %SETUP Sets up the linear system for the CMFD solve
            %   obj - cmfdClass to set up
            
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

