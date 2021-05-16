close all
tol = 1e-13;
lineSeg1 = [0,-0.5,1.1,0.5];
lineSeg2 = [0.1,0,1,0];

line_segment_intersect_indicator(lineSeg1,lineSeg2,tol)

lnwidth = 2;
plot(lineSeg1(1:2:3),lineSeg1(2:2:4),'r-','linewidth',lnwidth)
hold on
plot(lineSeg2(1:2:3),lineSeg2(2:2:4),'k-','linewidth',lnwidth)