%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan on 8/4/2013
%%% Last modified date: 6/28/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input module for finite elemnt 
%              and sensitivity analysis


function mesh = read_abaqus_mesh(inputfile)
    
%%% Introducing mesh parameters

%%%  NODE

% mesh.node.n_node: number of nodes in mesh
%mesh.node.node_n: node numebr
%mesh.node.coords: node coordinate
%mesh.node.Dirichlet.n_pre_temp: the number of prescribed nodal value
%mesh.node.Dirichlet.temp_node: the node number of prescribed nodal value
%mesh.node.Dirichlet.temp_value: the value of prescribed nodal value

%%% ELEMENT

%mesh.elem.n_elem: number of element in mesh
%mesh.elem.elem_n: element number
%mesh.elem.type: element type
%mesh.elem.elem_node: node connectivity 
%mesh.elem.heatSource: element heat source
%mesh.elem.material: element materail
%mesh.elem.Neumann.n_heatFlux: the number of element with Neumann BC
%mesh.elem.Neumann.heatFlux_elem: the element number of  element with Neumann BC
%mesh.elem.Neumann.heatFlux_value: the value of Neumann BC

%%% MATERIAL
%mesh.elem.material.name: material name
%mesh.elem.material.density: density
%mesh.elem.material.conductivity: conductivity
%mesh.elem.material.specificHeat: specific heat


%%%%BCs

%Dirichlet BCs
%mesh.BCs.Dirichlet(i_NSET).n_pre_temp: number of pescribed nodal
%                                                              temperature
%mesh.BCs.Dirichlet(i_NSET).temp_node: the node number belongs to BCs
%mesh.BCs.Dirichlet(i_NSET).temp_value: the prescribed nodal value

%Neumann BCs
%mesh.BCs.Neumann(i_Neumann_BCs).n_heatFlux: number of element subjected to
%                                                                          heat flux
%mesh.BCs.Neumann(i_Neumann_BCs).heatFlux_elem: the element number belongs
%                                                                               to BCs
%mesh.BCs.Neumann(i_Neumann_BCs).heatFlux_value: Heat flux on the Boundary

%%% Body Forces (Heat Sources)
%mesh.heatSource(i_heatSource).heatSource_elem: the element number
%mesh.heatSource(i_heatSource).heatSource_value: Heat source value



format long e;

% This function reads the input file.
fileID=fopen(inputfile);
%%%%%%%%%%%%%%%%%%%%%%%%%
%  Finite Element Input %
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% Read the input file: the first time
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initilize number of nodes, elements, node sets, and element sets
 mesh.node.n_node = 0; % number of nodes
 mesh.elem.n_elem = int32(0); % number of elements
 n_NSET = int32(0);  % Node set counter 
 n_ELSET = int32(0);  % Element set counter

 %for i=1:2000
while true % First while loop: 1st time reading
    line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
    if strcmp(line_string, '*END STEP')
        break;
    %end
    
%%%%%%%%%%%%%%%%%    
% Read number of node in the mesh
%%%%%%%%%%%%%%%%%

    elseif strcmp(line_string, '*NODE')
        while true % The Second while loop: 1st time reading
            line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
            if isempty(strfind(line_string, '**'))
                 mesh.node.n_node =  mesh.node.n_node + 1;
            else
                break;
            end
        end % end for the second while loop: 1st time reading

    %end of elseif strfind(line_string, '*NODE')

%%%%%%%%%%%%%%%%%%%
% Read number of element in the mesh
%%%%%%%%%%%%%%%%%%%

    elseif strfind(line_string, '*ELEMENT,')
        while true % The Third while loop: 1st time reading
            line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
            if isempty(strfind(line_string, '**'))
                 mesh.elem.n_elem =  mesh.elem.n_elem + 1;
            else
                break;
            end
        end % end for the third while loop: 1st time reading

    %end of elseif strfind(line_string, '*ELEMENT')

%%%%%%%%%%%%%%%%%%%    
% Read number of node sets in the mesh
%%%%%%%%%%%%%%%%%%%
    elseif strfind(line_string, '*NSET')
        n_NSET =  n_NSET + 1;
        
    %end of elseif strfind(line_string, '*NSET')    
   
%%%%%%%%%%%%%%%%%%%%%    
% Read number of Element sets in the mesh
%%%%%%%%%%%%%%%%%%%%%

    elseif strfind(line_string, '*ELSET')
        n_ELSET =  n_ELSET + 1;

    %end of elseif strfind(line_string, '*NSET')    

    
    end % end for if strcmp(line_string, '*END STEP')
    
end %end for while TRUE : the 1st time reading
 
%
%%%%%%%%%%%%%%%%%%%
% End of reading input file: the first time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% Read the input file: the second time
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frewind(fileID); % Back to the beginning of the input file

%%%%%%%%%%%%%%%%%%%%%%%%%
%initilize node-sets and element-sets name matrices
%%%%%%%%%%%%%%%%%%%%%%%%%
if (n_NSET > 0)
    node_set(n_NSET).name = [];
end
if (n_ELSET > 0)
    elem_set(n_ELSET).name = [];
end
i_Neumann_BCs = int32(0); % Neumann BCs counter
i_heatSource = int32(0); % Heat Source counter

%%%%%%%%%%%%%
%%%Some needed counter
%%%%%%%%%%%%%

i_NSET = int32(0);  % Node-set counter
i_ELSET = int32(0); %Element-set counter


while true % First while loop: the 2nd time reading
    line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
    if strcmp(line_string, '*END STEP')
        break;
    %end
        
%%%%%%%%%%%%%
%Read coordinates of nodes
%%%%%%%%%%%%%

    elseif strcmp(line_string, '*NODE')
        for i = 1: mesh.node.n_node
            line = fscanf(fileID, '%i, %f, %f\n', 3);
            mesh.node.node_n(i) = int32(line(1));
            mesh.node.coords(mesh.node.node_n(i), 1:2) = line(2:3);
        end
    %end of elseif strfind(line_string, 'NODE')
    
%%%%%%%%%%%%%
%Read elements information
%%%%%%%%%%%%%
   
    elseif strfind(line_string, '*ELEMENT,')
        %read element type and element set
        %{
        split_string = regexp(line_string, ',', 'split');
        split_string2 = regexp(split_string(2), '=', 'split');
        mesh.elem.type = split_string2{1, 1}(1, 2);
        split_string2 = regexp(split_string(3), '=', 'split');
        mesh.material_properties_name = split_string2{1, 1}(1, 2);
        %}
        %read element conectivity
        for i = 1: mesh.elem.n_elem
            line = fscanf(fileID, '%i, %i, %i, %i\n', 4);
            mesh.elem.elem_n(i) = int32(line(1));
            mesh.elem.elem_node(mesh.elem.elem_n(i),1:3) = line(2:4);
        end
    %end of elseif strfind(line_string, '*ELEMENT')
    
    elseif strfind(line_string, '*SHELL SECTION')
         split_string = regexp(line_string, ',', 'split');
         split_string2 = regexp(split_string(3), '=', 'split');
         mesh.material.name = split_string2{1, 1}(1, 2);
         
         %read shell thickness and number of devision
         %{
         line = fscanf(fileID, '%f, %f\n', 2);
         mesh.geometry.shell_thickness = line (1);
         mesh.geometry.shell_devision = line(2);
         %}
     %end of elseif strfind(line_string, '*SHELL SECTION')

%%%%%%%%%%%%%
% read material properties
%%%%%%%%%%%%%

    elseif strfind(line_string, '*MATERIAL')
          split_string1 = regexp(line_string, ',', 'split');
         material_name1 = regexp(split_string1(2), '=', 'split');
         if strcmp(mesh.material.name, material_name1{1, 1}(1, 2))
             while true % second while loop
                 line_string1 = fgetl(fileID);
                 if strfind(line_string1, '** Step  1')
                     break;
                 end
                 
                 if strfind(line_string1, 'DENSITY')
                     line = fscanf(fileID, '%f,\n', 1);
                     mesh.material.density = line (1);
                 %end
                 
                 elseif strfind(line_string1, 'CONDUCTIVITY')
                      split_string = regexp(line_string1, ',', 'split');
                      type = regexp(split_string(2), '=', 'split');
                      if strcmp('ISO', type{1, 1}(1, 2))
                          line = fscanf(fileID, '%f,\n', 1);
                          mesh.material.conductivity = line (1);
                      end
                  %end  % end for if strfind(line_string1, 'CONDUCTIVITY')
                  
                 elseif strfind(line_string1, 'SPECIFIC HEAT')
                      line = fscanf(fileID, '%f,\n', 1);
                      mesh.material.specificHeat = line (1);
                 end   %end for if strfind(line_string1, 'DENSITY')
                 
             end % end for second while loop
             
         end % end for if strcmp(elem_material_name, material_name{1, 1}(1, 2));
      
     % end of elseif strfind(line_string, '*MATERIAL')
     
%%%%%%%%%%%%%%%%%%%%%%%%%%
% read B.Cs. and body forces node-set and element-set
%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    elseif strfind(line_string, '*NSET')
        n_node_in_NSET = 0; % initialize the number of node in the current node_set to zero
        i_NSET = i_NSET +1; % adding one to element_set counter
        %read Node Set (Dirichlet B.C.)
        split_string = regexp(line_string, ',', 'split');
        split_string2 = regexp(split_string(2), '=', 'split');
        node_set(i_NSET).name = split_string2{1, 1}(1, 2);
        SIZE = size(split_string);
        if SIZE(2) > 2 
            if strcmp(split_string(3), ' GENERATE')
                while true % The  second while loop: the 2nd time reading
                    position_1 = ftell(fileID); % Read current position of file
                    line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                    position_2 = ftell(fileID); % Read current position of file
                    % Back to previuos line and do nothing
                    fseek(fileID, (position_1-position_2),'cof'); % Back to previuos line
                    %position_3 = ftell(fileID); % Read current position of file
                    if strfind(line_string, '*')
                        break;
                    else
                        line = fscanf(fileID, '%i, %i, %i\n',3);
                        s_n = int32(line(1)); % Starting node for BC1 with prescribe temperature
                        e_n = int32(line(2)); % Ending node for BC1 with prescribe temperature
                        step = int32(line(3)); % Step from starting node to ending node
                        node_set(i_NSET).number = (e_n-s_n)/step + 1;
                        node_set(i_NSET).nodes = int32(zeros(node_set(i_NSET).number, 1));
                        for i =1:node_set(i_NSET).number
                            node_set(i_NSET).nodes( i) = s_n + (i - 1) * step;
                        end
                    end
                end % end for the  second while loop: the 2nd time reading
            end % end for strcmp(split_string(3), 'GENERATE')

        elseif  SIZE(2) == 2
            position_1 = ftell(fileID); % Read current position of file
            NSET_line_number = 0; % initialize node-set-line counter
            while true %The third while loop: the 2nd time reading
                line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                position_2 = ftell(fileID); % Read current position of file
                if strfind(line_string, '*')
                    fseek(fileID, (position_1-position_2),'cof'); % Back to first line in elem_set
                    node_set(i_NSET).number = n_node_in_NSET; % asignning the number of node in the current node-set
                    node_set(i_NSET).nodes = int32(zeros(node_set(i_NSET).number, 1)); % initialize the node matrix in the current node-set 
                    start = 1; % initialize set-node counter 
                    for line_number = 1:NSET_line_number
                        line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                        split_string = regexp(line_string, ',', 'split');
                        SIZE3 = size(split_string);
                        if isempty(char(split_string(1, end)))
                            SIZE3(2) = SIZE3(2) - 1;
                        end
                        for i = start:(SIZE3(2) + start - 1)
                            node_set(i_NSET).nodes(i) = str2double(split_string(i - start + 1)); % Store the node number in the set
                        end
                        start = i+1;
                        
                    end % end of "for line_number = 1:NSET_line_number"
                    break;
                    
                else
                    NSET_line_number = NSET_line_number +1;
                    split_string = regexp(line_string, ',', 'split');
                    SIZE2 = size(split_string);
                    %AAA = char(split_string(1, end));
                    if isempty(char(split_string(1, end)))
                        SIZE2(2) = SIZE2(2) - 1;
                    end
                    n_node_in_NSET = n_node_in_NSET + SIZE2(2);
                 
                end % end for if strfind(line_string, '*')
            end %The third while loop: the 2nd time reading
        end %end for if SIZE(2) > 2
        
        %end of elseif strfind(line_string, '*NSET')
        
    elseif strfind(line_string, '*ELSET')
        n_elem_in_ELSET = 0; % initialize the number of element in the current element_set to zero
        i_ELSET = i_ELSET +1; % adding one to element_set counter
        %read Node Set (Dirichlet B.C.)
        split_string = regexp(line_string, ',', 'split');
        split_string2 = regexp(split_string(2), '=', 'split');
        elem_set(i_ELSET).name = split_string2{1, 1}(1, 2);
        SIZE = size(split_string);
        if SIZE(2) > 2 
            if strcmp(split_string(3), ' GENERATE')
                while true % The forth while loop: the 2nd time reading
                    position_1 = ftell(fileID); % Read current position of file
                    line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                    position_2 = ftell(fileID); % Read current position of file
                    % Back to previuos line and do nothing
                    fseek(fileID, (position_1-position_2),'cof'); % Back to previuos line
                    %position_3 = ftell(fileID); % Read current position of file
                    if strfind(line_string, '*')
                        break;
                    else
                        line = fscanf(fileID, '%i, %i, %i\n',3);
                        s_n = line(1); % Starting node for BC1 with prescribe temperature
                        e_n = line(2); % Ending node for BC1 with prescribe temperature
                        step = line(3); % Step from starting node to ending node
                        elem_set(i_ELSET).number = (e_n-s_n)/step + 1; % asignning the number of element in the current element-set
                        elem_set(i_ELSET).elem = int32(zeros(elem_set(i_ELSET).number, 1)); % initialize the element matrix in the current element-set 
                        for i =1:elem_set(i_ELSET).number
                            elem_set(i_ELSET).elem( i) = s_n + (i - 1) * step;
                        end
                    end
                end % end for the forth while loop: the 2nd time reading
            end % end for strcmp(split_string(3), 'GENERATE')
            
        elseif  SIZE(2) == 2
            position_1 = ftell(fileID); % Read current position of file
            ELSET_line_number = 0; % initialize element-set-line counter
            while true %The fifth while loop: the 2nd time reading
                line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                position_2 = ftell(fileID); % Read current position of file
                if strfind(line_string, '*')
                    fseek(fileID, (position_1-position_2),'cof'); % Back to first line in elem_set
                    elem_set(i_ELSET).number = n_elem_in_ELSET; % asignning the number of element in the current element-set
                    elem_set(i_ELSET).elem = zeros(elem_set(i_ELSET).number, 1); % initialize the element matrix in the current element-set 
                    start = 1; % initialize set-element counter 
                    for line_number = 1:ELSET_line_number
                        line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
                        split_string = regexp(line_string, ',', 'split');
                        SIZE3 = size(split_string);
                        if isempty(char(split_string(1, end)))
                            SIZE3(2) = SIZE3(2) - 1;
                        end
                        for i = start:(SIZE3(2) + start - 1)
                            elem_set(i_ELSET).elem(i) = str2double(split_string(i - start + 1)); % Store the element number in the set
                        end
                        start = i+1;
                        
                    end % end of "for line_number = 1:ELSET_line_number"
                    break;
                    
                else
                    ELSET_line_number = ELSET_line_number +1;
                    split_string = regexp(line_string, ',', 'split');
                    SIZE2 = size(split_string);
                    %AAA = char(split_string(1, end));
                    if isempty(char(split_string(1, end)))
                        SIZE2(2) = SIZE2(2) - 1;
                    end
                    n_elem_in_ELSET = n_elem_in_ELSET + SIZE2(2);
                 
                end % end for if strfind(line_string, '*')
            end %The fifth while loop: the 2nd time reading
        end %end for if SIZE(2) > 2
        
        %end of elseif strfind(line_string, '*ELSET')

%%%%%%%%%%%%%%%%%%
% Store B.Cs. and body forces values
%%%%%%%%%%%%%%%%%%
    
    % Dirichlet BCs
    elseif strfind(line_string, '*BOUNDARY')
        line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
        split_string = regexp(line_string, ',', 'split');
        split_string(1)=strtrim(split_string(1));
        for i_NSET = 1 : n_NSET
            if strcmp(split_string(1), node_set(i_NSET).name)
                mesh.BCs.Dirichlet(i_NSET).temp_node = node_set(i_NSET).nodes(:);
                mesh.BCs.Dirichlet(i_NSET).n_pre_temp = node_set(i_NSET).number;
                for i_node_in_NSET = 1:node_set(i_NSET).number
                    mesh.BCs.Dirichlet(i_NSET).temp_value(i_node_in_NSET, 1) = str2double(split_string(end));
                end
            end
        end
   %end of elseif strfind(line_string, '*BOUNDARY')

   % Neumann BCs and Body Forces
    elseif strfind(line_string, '*DFLUX')
        line_string = fgetl(fileID);   % line_string is a dummy variable that stores string associated with each line.
        split_string = regexp(line_string, ',', 'split');
        split_string(1)=strtrim(split_string(1));
        split_string(2)=strtrim(split_string(2));
        for i_ELSET = 1:n_ELSET
            if strcmp(strtrim(split_string(1)), elem_set(i_ELSET).name)
                %Neumann BCs
                if strcmp(split_string(2), 'S1') || strcmp(split_string(2), 'S2') || strcmp(split_string(2), 'S3')
                    i_Neumann_BCs = i_Neumann_BCs +1;
                    mesh.BCs.Neumann(i_Neumann_BCs).heatFlux_surface = split_string(2);
                    mesh.BCs.Neumann(i_Neumann_BCs).heatFlux_elem = elem_set(i_ELSET).elem(:);
                    mesh.BCs.Neumann(i_Neumann_BCs).n_heatFlux = elem_set(i_ELSET).number;
                    for i_elem_in_ELSET = 1 : elem_set(i_ELSET).number
                        mesh.BCs.Neumann(i_Neumann_BCs).heatFlux_value(i_elem_in_ELSET, 1) = str2double(split_string(end));
                    end
                    
                %Body Forces
                elseif strcmp(split_string(2), 'BF')
                    i_heatSource = i_heatSource +1;
                    %mesh.elem.heatSource.heatSource_elem = elem_set(i_ELSET).elem(:);
                    heatSource_elem = elem_set(i_ELSET).elem(:);
                    for i_elem_in_ELSET = 1 : elem_set(i_ELSET).number
                        %mesh.elem.heatSource.heatSource_value(i_elem_in_ELSET, 1) = str2double(split_string(end));
                        heatSource_value(i_elem_in_ELSET, 1) = str2double(split_string(end));
                    end
                    %mesh.heatSource(i_heatSource).heatSource_value = str2double(split_string(end));
                    %mesh.heatSource(i_heatSource).heatSource_elem = elem_set(i_ELSET).elem(:);
                end
            end
        end
    %end of elseif strfind(line_string, '*DFLUX')
   
    end % end for if strcmp(line_string, '*END STEP')
             
end %end for the first while TRUE: the 2nd time reading


%%%Assign material properties to each element
%{
mesh.elem.material.density = zeros(mesh.elem.n_elem, 1);
mesh.elem.material.conductivity = zeros(mesh.elem.n_elem, 1);
mesh.elem.material.specificHeat = zeros(mesh.elem.n_elem, 1);

for i = 1: mesh.elem.n_elem
            mesh.elem.material.density(i) = mesh.material.density;
            mesh.elem.material.conductivity(i) = mesh.material.conductivity;
            mesh.elem.material.specificHeat(i) = mesh.material.specificHeat;
end
%}
%%%%%%%%%%
warning('only one kind of material is allowed');
mesh.elem.material(1: mesh.elem.n_elem) = int32(1);
%%%%%%%%%
%%% store all the Dirichlet BCs in a structure
mesh.node.Dirichlet.n_pre_temp = int32(0); % initialize number of prescribed nodal value to zero
mesh.node.Dirichlet.temp_node = [];
mesh.node.Dirichlet.temp_value = [];

mesh.elem.Neumann.n_heatFlux = int32(0);
mesh.elem.Neumann.heatFlux_elem = [];
mesh.elem.Neumann.heatFlux_value = [];
mesh.elem.Neumann.heatFlux_surface = [];

if(isfield(mesh,'BCs'))
    %%% store all the Dirichlet BCs in a structure
    if(isfield(mesh.BCs,'Dirichlet'))
        n_Dirichlet = size(mesh.BCs.Dirichlet, 2); % The number of Dirichlet BCs
    end
   
    for i_Dirichlet = 1 : n_Dirichlet
        mesh.node.Dirichlet.n_pre_temp = mesh.node.Dirichlet.n_pre_temp + mesh.BCs.Dirichlet(i_Dirichlet).n_pre_temp;
    end

    mesh.node.Dirichlet.temp_node = int32(zeros(mesh.node.Dirichlet.n_pre_temp, 1)); % initialize the matrix which stores the node number of prescribed nodal value 
    mesh.node.Dirichlet.temp_value = zeros(mesh.node.Dirichlet.n_pre_temp, 1); % initialize the matrix which stores the prescribed nodal value

    k = 0;
    for i_Dirichlet = 1 : n_Dirichlet
        j = k + 1;
        k = k + mesh.BCs.Dirichlet(i_Dirichlet).n_pre_temp;
        mesh.node.Dirichlet.temp_node(j:k)  = mesh.BCs.Dirichlet(i_Dirichlet).temp_node(:);
        mesh.node.Dirichlet.temp_value(j:k)  = mesh.BCs.Dirichlet(i_Dirichlet).temp_value(:);
    end


    %%% store all the Neumann BCs in a structure
    if(isfield(mesh.BCs,'Neumann'))
        n_Neumann = size(mesh.BCs.Neumann, 2); % The number of Neumann BCs
        mesh.elem.Neumann.n_heatFlux = int32(0); % initialize number of element with Neumann BC to zero

        for i_Neumann = 1 : n_Neumann
            mesh.elem.Neumann.n_heatFlux = mesh.elem.Neumann.n_heatFlux + mesh.BCs.Neumann(i_Neumann).n_heatFlux;
        end
        
        mesh.elem.Neumann.heatFlux_elem = int32(zeros(mesh.elem.Neumann.n_heatFlux, 1)); % initialize the matrix which stores the element number of element with Neumann BC 
        mesh.elem.Neumann.heatFlux_value = zeros(mesh.elem.Neumann.n_heatFlux, 1); % initialize the matrix which stores the Neumann BC value
        mesh.elem.Neumann.heatFlux_surface = repmat(cellstr(''), mesh.elem.Neumann.n_heatFlux, 1);
        
        k = 0;
        for i_Neumann = 1 : n_Neumann
            j = k + 1;
            k = k + mesh.BCs.Neumann(i_Neumann).n_heatFlux;
            mesh.elem.Neumann.heatFlux_surface(j:k) = mesh.BCs.Neumann(i_Neumann).heatFlux_surface;
            mesh.elem.Neumann.heatFlux_elem(j:k)  = mesh.BCs.Neumann(i_Neumann).heatFlux_elem(:);
            mesh.elem.Neumann.heatFlux_value(j:k)  = mesh.BCs.Neumann(i_Neumann).heatFlux_value(:);
        end
    %else
    %    mesh.BCs.Neumann=[];
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh.elem.heatSource = zeros(mesh.elem.n_elem,1);
mesh.elem.heatSource(heatSource_elem) = heatSource_value;
for i_ELSET = 1:n_ELSET
    sInd = cell2mat(regexp(elem_set(i_ELSET).name,'REGION','once'));
    if(~isempty(sInd))         
        regionNum = str2double(elem_set(i_ELSET).name{1}(7));
        mesh.elem.region(elem_set(i_ELSET).elem) = regionNum; 
    end
end

mesh.elem = rmfield(mesh.elem,'elem_n');
mesh.node = rmfield(mesh.node,'node_n');
mesh.boundary.xi = min(mesh.node.coords(:,1));
mesh.boundary.xf = max(mesh.node.coords(:,1));
mesh.boundary.yi = min(mesh.node.coords(:,2));
mesh.boundary.yf = max(mesh.node.coords(:,2));
mesh = rmfield(mesh,'BCs');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fileID);

return
