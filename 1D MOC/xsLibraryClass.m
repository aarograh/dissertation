classdef xsLibraryClass < handle
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
        function obj = xsLibraryClass( input )
            %XSLIBRARYCLASS Constructor for xsLibraryClass
            %   input - The inputClass container from which to initialize
            
            obj.openfile(input.xsfilename);
            
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
                while isempty(k) || nextline(1) == '!'
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
                % Calculate total cross-section and transport cross-section
                obj.xsSets(i).calcTXS( input.scattype );
            end
            
            obj.closefile( );
            
            for j=1:input.nmixtures
                id = input.mixtures(j,1);
                name = sprintf('Mixture ');
                for k=2:length(input.mixtures(j,:))
                    if input.mixtures(j,k) == 0
                        break
                    else
                        if k > 2
                            name = sprintf('%s;',name);
                        end
                        name = sprintf('%s %0.3f%% %s',name,input.mixvols(j,k-1),...
                            obj.xsSets(input.mixtures(j,k)).name);
                    end
                end
                obj.xsSets(id) = xsClass(name);
                obj.xsSets(id).nsubxs = 0;
                obj.xsSets(id).subxs = xsClass();
                for k=2:length(input.mixtures(j,:))
                    if input.mixtures(j,k) == 0
                        break
                    else
                        obj.xsSets(id).nsubxs = obj.xsSets(id).nsubxs + 1;
                        obj.xsSets(id).subxs(k-1) = obj.xsSets(input.mixtures(j,k));
                    end
                end
                obj.xsSets(id).subfracs = input.mixvols(j,:);
                obj.mix( id );
            end
        end
        
        function obj = mix( obj, imix )
            %MIX Mixes sub-cross-sections
            %   obj     - the xsLibraryClass object
            %   imix    - the index of the mixture being set up
            
            % Initialize cross-sections and get scattering order
            obj.xsSets(imix).scatOrder = max(obj.xsSets(imix).subxs.scatOrder);
            obj.xsSets(imix).total(1:obj.ngroups) = 0;
            obj.xsSets(imix).transport(1:obj.ngroups) = 0;
            obj.xsSets(imix).absorption(1:obj.ngroups) = 0;
            obj.xsSets(imix).nufission(1:obj.ngroups) = 0;
            obj.xsSets(imix).fission(1:obj.ngroups) = 0;
            obj.xsSets(imix).chi(1:obj.ngroups) = 0;
            obj.xsSets(imix).scatter(1:obj.ngroups,1:obj.ngroups, ...
                obj.xsSets(imix).scatOrder+1) = 0;
            for i=1:obj.xsSets(imix).nsubxs
                
                obj.xsSets(imix).total = obj.xsSets(imix).total + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).total;
                obj.xsSets(imix).transport = obj.xsSets(imix).transport + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).transport;
                obj.xsSets(imix).absorption = obj.xsSets(imix).absorption + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).absorption;
                obj.xsSets(imix).nufission = obj.xsSets(imix).nufission + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).nufission;
                obj.xsSets(imix).fission = obj.xsSets(imix).fission + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).fission;
                obj.xsSets(imix).chi = obj.xsSets(imix).chi + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).chi;
                obj.xsSets(imix).scatter = obj.xsSets(imix).scatter + ...
                    obj.xsSets(imix).subfracs(i)*obj.xsSets(imix).subxs(i).scatter;
            end
        end
        
        function obj = openfile( obj, filename )
            %OPENFILE Opens the XS Library file
            %   obj      - The xsLibraryClass object
            %   filename - The XS Library filename
            
            obj.fileid = fopen(filename);
            
        end
        
        function obj = closefile( obj )
            %CLOSEFILE Closes the XS Library file
            %   obj      - The xsLibraryClass object
            
            fclose(obj.fileid);
            obj.fileid = 0;
            
        end
        
        function [ nextline, EOF ] = getLine( obj )
            %GETLINE Gets the next line from the XS Library file
            %   obj      - The xsLibraryClass object
            
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

