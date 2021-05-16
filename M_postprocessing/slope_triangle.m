% this function generates a slope triangle in the log log plot given the 
% left-most point of the flat edge of the triangle (x1,y1), the slope of 
% the slanted edge of the triangle and the length of the flat edge
function [xpt,ypt]=slope_triangle(xdata,ydata,slope)
[xdata,ix]=sort(xdata);
ydata = ydata(ix);
if(slope<0)
    [~,ind] = max(xdata);
    ind = ind-1;    
else
    [~,ind] = min(xdata);
end
h = xdata(ind+1)-xdata(ind);
x1 = xdata(ind);
yoffset = ydata(ind+1)-ydata(ind);
y1 = ydata(ind)+yoffset;
xpt(1:4)=0;
ypt(1:4)=0;
xpt(1) = x1;
ypt(1) = y1;
xpt(2) = x1+h;
ypt(2) = y1*(xpt(2)/xpt(1))^slope;
xpt(3) = x1;
ypt(3) = ypt(2);
xpt(4)=xpt(1);
ypt(4)=ypt(1);

end