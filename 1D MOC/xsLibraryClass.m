classdef xsLibraryClass
    %XSLIBRARYCLASS Stores cross-section data and reads XS library file
    %   Cross-section libraries must be in the same format as the
    %   MPACT user libraries.  One cross-section set will be stored for
    %   each material type.
    
    properties
        fileid
        name
        ngroups
        nsets
        groupBounds
        xsSets
    end
    
    methods
        function obj = xsLibraryClass( filename )
            obj = obj.openfile(filename);
            
            % Get library name
            nextline = obj.getLine();
            obj.name = strrep(nextline,sprintf('\n'),'');
            
            % Get number of energy groups and materials
            nextline = obj.getLine();
            [obj.ngroups, obj.nsets] = strread(nextline);
            
            % Get energy group boundaries
            nextline = obj.getLine();
            obj.groupBounds = strread(nextline);
            
            for i=1:obj.nsets
                % Loop until we actually hit the XS set
                nextline = obj.getLine();
                k = strfind(nextline,'XSMACRO');
                while isempty(k)
                    nextline = obj.getLine();
                    k = strfind(nextline,'XSMACRO');
                end
                
                % Read name and number of scattering moments
                tmpstr = strsplit(nextline);
                % Split the stupid cell array object into string and integer
                setname = strtrim(cell2mat(tmpstr(2)));
                setorder = strread(strtrim(cell2mat(tmpstr(3))));
                % This branching is needed to prevent matlab from assuming xsSets is a double array
                if i == 1
                    obj.xsSets = xsClass(setname, setorder);
                else
                    obj.xsSets(i) = xsClass(setname, setorder);
                end
                
                % Read in absorption, nufission, fission, and chi
                for k=1:obj.ngroups
                    nextline = obj.getLine();
                    [obj.xsSets(i).absorption(k), obj.xsSets(i).nufission(k), ...
                        obj.xsSets(i).fission(k), obj.xsSets(i).chi(k)] = strread(nextline);
                end
                % Loop over scattering orders
                for j=1:setorder+1
                    % Read in scattering matrix for order j-1
                    for k=1:obj.ngroups
                        nextline = obj.getLine();
                        obj.xsSets(i).scatter(k,1:obj.ngroups,j) = strread(nextline);
                    end
                end
            end
            
            obj = obj.closefile();
        end
        
        function obj = openfile( obj, filename )
            obj.fileid = fopen(filename);
        end
        
        function obj = closefile( obj )
            fclose(obj.fileid);
            obj.fileid = 0;
        end
        
        function [ nextline, EOF ] = getLine( obj )
            nextline = fgets(obj.fileid);
            if ischar(nextline)
                nextline = strrep(nextline,sprintf('\n'),'');
                EOF = 0;
            else
                EOF = 1;
            end
        end
    end
    
end

