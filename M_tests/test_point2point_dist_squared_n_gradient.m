path(path,'../M_geom_toolbox')
rng('shuffle')
pts1 = rand(4,2);
pts2 = rand(4,2);

[distSq,distSqGrad] = point2point_dist_squared_n_gradient(pts1,pts2);

del = 1e-7;
distSqFD = nan(size(distSqGrad));

newPts1 = pts1;
newPts1(:,1) = newPts1(:,1) - del; 
newDistSq1 = point2point_dist_squared_n_gradient(newPts1,pts2);
newPts1 = pts1;
newPts1(:,1) = newPts1(:,1) + del;
newDistSq2 = point2point_dist_squared_n_gradient(newPts1,pts2);
distSqFD(:,1) = 0.5*(newDistSq2 - newDistSq1)/del; 

newPts1 = pts1;
newPts1(:,2) = newPts1(:,2) - del; 
newDistSq1 = point2point_dist_squared_n_gradient(newPts1,pts2);
newPts1 = pts1;
newPts1(:,2) = newPts1(:,2) + del;
newDistSq2 = point2point_dist_squared_n_gradient(newPts1,pts2);
distSqFD(:,2) = 0.5*(newDistSq2 - newDistSq1)/del; 


newPts2 = pts2;
newPts2(:,1) = newPts2(:,1) - del; 
newDistSq1 = point2point_dist_squared_n_gradient(pts1,newPts2);
newPts2 = pts2;
newPts2(:,1) = newPts2(:,1) + del;
newDistSq2 = point2point_dist_squared_n_gradient(pts1,newPts2);
distSqFD(:,3) = 0.5*(newDistSq2 - newDistSq1)/del; 

newPts2 = pts2;
newPts2(:,2) = newPts2(:,2) - del; 
newDistSq1 = point2point_dist_squared_n_gradient(pts1,newPts2);
newPts2 = pts2;
newPts2(:,2) = newPts2(:,2) + del;
newDistSq2 = point2point_dist_squared_n_gradient(pts1,newPts2);
distSqFD(:,4) = 0.5*(newDistSq2 - newDistSq1)/del; 


relDiff = abs(1 - distSqFD./distSqGrad);

fprintf('max relative difference = %g \n', max(relDiff(:)));
