%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/6/2013
%%% Last modified date: 8/6/2013
%%% Copyright 2013 University of Illinois 
function edge_length=find_edge_length(edgeNode,nodeCoord)
    nEdge=size(edgeNode,1);
    edge_length=zeros(nEdge,1);
    for i=1:nEdge
        r1=nodeCoord(edgeNode(i,1),:);
        r2=nodeCoord(edgeNode(i,2),:);
        edge_length(i)=norm(r2-r1,2);
    end
end