%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/17/2013
%%% Last modified date: 9/17/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function plot the convexhulls of nurbs curve
function plot_nurbs_convexhull(itrface)
subd=50;
figure
hold on
for i=1:itrface.nNurbs
    nurbs_plot(itrface.nurbs(i),subd);
    for k=1:itrface.nurbs(i).nConvHull
        X=[itrface.nurbs(i).hull{k}(1,:),itrface.nurbs(i).hull{k}(1,1)];
        Y=[itrface.nurbs(i).hull{k}(2,:),itrface.nurbs(i).hull{k}(2,1)];
        plot(X,Y,'k--','linewidth',1)
        text(mean(X),mean(Y),[num2str(i),',',num2str(k)])
    end
end
end