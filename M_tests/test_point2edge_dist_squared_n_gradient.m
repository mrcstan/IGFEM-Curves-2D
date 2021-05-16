path(path,'../M_geom_toolbox')
rng('shuffle')
pts = rand(4,2);
edges = rand(4,4);

[distSq,distSqGrad] = point2edge_dist_squared_n_gradient(pts,edges);

del = 1e-7;
distSqFD = nan(size(distSqGrad));

for i = 1:2
    newPts = pts;
    newPts(:,i) = newPts(:,i) - del; 
    newDistSq1 = point2point_dist_squared_n_gradient(newPts,edges);
    newPts = pts;
    newPts(:,i) = newPts(:,i) + del;
    newDistSq2 = point2point_dist_squared_n_gradient(newPts,edges);
    distSqFD(:,i) = 0.5*(newDistSq2 - newDistSq1)/del; 
end

for i = 1:4
    newEdges = edges;
    newEdges(:,i) = newEdges(:,i) - del; 
    newDistSq1 = point2point_dist_squared_n_gradient(pts,newEdges);
    newEdges = edges;
    newEdges(:,i) = newEdges(:,i) + del;
    newDistSq2 = point2point_dist_squared_n_gradient(pts,newEdges);
    distSqFD(:,i+2) = 0.5*(newDistSq2 - newDistSq1)/del; 
end

abdDiff = abs(distSqFD - distSqGrad);
relDiff = abs(1 - distSqFD./distSqGrad);

fprintf('max absolute difference = %g \n', max(absDiff(:)));
fprintf('max relative difference = %g \n', max(relDiff(:)));
