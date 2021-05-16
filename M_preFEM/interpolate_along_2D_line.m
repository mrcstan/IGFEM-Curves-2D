%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/1/2013
%%% Last modified date: 9/1/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function interpolates along a line in 2D when the value at the two
% end points of the line is specified.
% ASSUMPTION: the points (xI,yI) are on the line 
function val=interpolate_along_2D_line(x,y,z,xI,yI)
nx=length(x);
ny=length(y);
nz=length(z);
if(nx~=2 || ny~=2 || nz~=2)
    error('interpolate_along_2D_line: requires length(x)=length(y)=length(z)=2')
end
nxI=length(xI);
nyI=length(yI);
if(nxI~=nyI)
    error('interpolate_along_2D_line: requires length(xI)=length(yI)')
end
val=zeros(nxI,1);
totDist=norm([x(2)-x(1),y(2)-y(1)]);
delta=z(2)-z(1);
for i=1:nxI
    dist=norm([xI(i)-x(1),yI(i)-y(1)]);
    val(i)=z(1)+dist/totDist*delta;
end
end