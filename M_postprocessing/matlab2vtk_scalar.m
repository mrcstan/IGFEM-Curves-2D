%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/4/2014
%%% Last modified date: 2/26/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function converts a matlab output to a vtk unstructured grid consisting
% of triangular elements data format
%Output file name
function matlab2vtk_scalar(outfile,scalarname,elem,node,UUR2)

    fid = fopen(outfile, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK from Matlab\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');   
    fprintf(fid, ['POINTS ' num2str(node.nOriginalEnrichNode) ' float\n']);
     
    fprintf(fid, '%f %f %f \n',[node.coords(1:node.nOriginalEnrichNode,:),zeros(node.nOriginalEnrichNode,1)]');
    
    %{
    if(isfield(elem.parent,'nChild'))
        nChilds = sum(cat(1,elem.parent(:).nChild));
    else
        nChilds = 0;
    end
    
    nRegular = sum(cat(2,elem.parent(:).type) < 2);
    totElems = nChilds+nRegular;
    %}
    totElems = 0;
    listSize = 0;
    for i = 1:elem.n_elem
        if(elem.parent(i).type > 1)
            for j = 1:numel(elem.parent(i).child)
                if(elem.parent(i).child(j).isTriangle)
                    listSize = listSize + 4;
                else
                    listSize = listSize + 5;
                end
                totElems = totElems + 1;
            end
        else
            listSize = listSize + 4;
            totElems = totElems + 1;
        end
        
    end
    % listSize = listSize+nRegular*4;
    fprintf(fid, ['\nCELLS ', num2str(totElems), ' ', num2str(listSize),'\n']);
    
    cellTypes = 5*ones(1,totElems);
    count = 0;
    for i = 1:elem.n_elem
        if(elem.parent(i).type > 1)
            for j = 1:numel(elem.parent(i).child)              
                count = count+1;
                if(elem.parent(i).child(j).isTriangle)
                    fprintf(fid, '%i %i %i %i \n',[3,elem.parent(i).child(j).nodes(1:3)-1]);
                else
                    fprintf(fid, '%i %i %i %i %i\n',[4,elem.parent(i).child(j).nodes(1:4)-1]); 
                    cellTypes(count) = 9;                   
                end
                
            end
        else
            count = count+1;
            fprintf(fid, '%i %i %i %i \n',[3,elem.elem_node(i,:)-1]);

        end
    end
  
    fprintf(fid, ['\nCELL_TYPES ',num2str(totElems),'\n']);
    fprintf(fid, '%i\n',cellTypes);
    fprintf(fid, ['\nPOINT_DATA ',num2str(node.nOriginalEnrichNode),'\n']);
    strg = ['SCALARS ',scalarname,' float 1\n'];
    fprintf(fid, strg);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f \n', UUR2(1:node.nOriginalEnrichNode));
    
    fclose(fid);
end