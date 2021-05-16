clear all
close all
path(path,'../M_geom_toolbox')
coordBounds = [0,1.0;0,1.0];
lengthBounds = [0.05,0.1;0.45,0.5;0.05,0.1;0.45,0.5];
minAngle = 5*pi/180;
%fixedPoints = [0.1,0.7;0.1,0.2];
%fixedPoints = [0.1,0.7];
fixedPoints = [];
maxTrial = 100;
rep = 1;
attempts = nan(rep,1);
for i = 1:rep
    [quad,~,attempts(i)] = generate_random_quad(coordBounds, ...
                                                lengthBounds, ...
                                                minAngle, ...
                                                fixedPoints, ...
                                                maxTrial, ...
                                                true);
end
fprintf('average number of attempts = %i \n',mean(attempts))
fprintf('max number of attempts = %i \n',max(attempts))
plot([quad(:,1);quad(1,1)],[quad(:,2);quad(1,2)],'r^-',...
     'linewidth',3,'markerfacecolor','r','markersize',12)
 
for i = 1:4
    text(quad(i,1),quad(i,2),num2str(i),'fontsize',20)
end