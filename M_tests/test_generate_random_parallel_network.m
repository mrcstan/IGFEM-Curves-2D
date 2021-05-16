clear all
close all
path(path,'../M_geom_toolbox')

nBranches = 5;
coordBounds = [0,0.15;0,0.2];
length1 = 0.08/nBranches;
length2 = 0.12/nBranches;
lengthBounds = [0.1,0.14;
                length1,length2;
                0.1,0.14;
                length1,length2];
minAngle = 5*pi/180;
maxTrial = 500;

iniVertex = [0,0.19];
finalVertex = [0.15,0.01];
[vertices,connectivity] = generate_random_parallel_network(iniVertex, ...
                                                           finalVertex, ...
                                                           nBranches, ...
                                                           coordBounds, ...
                                                           lengthBounds, ...
                                                           minAngle, ...
                                                           maxTrial, ...
                                                           true);
figure
hold on
nConnectivity = size(connectivity,1); 
for i = 1:nConnectivity
    if (i == 1 || i == nConnectivity)
        XX = vertices(connectivity(i,1:2),:);
    else
        XX = [vertices(connectivity(i,:),:);vertices(connectivity(i,1),:)];
    end
    plot(XX(:,1),XX(:,2),'-r','linewidth',3)
end

ftsz = 24;
axis image
xlabel('x','fontsize',ftsz)
xlabel('y','fontsize',ftsz)
set(gca,'fontsize',ftsz)