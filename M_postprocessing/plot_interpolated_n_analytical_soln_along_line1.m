%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/14/2013
%%% Last modified date: 1/22/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function plot the solution along a straight line specified by the
% array line = [xo,yo,vx,vy] where (xo,yo) is the starting point
% and (vx,vy) is the unit vector parallel to the line
function varargout = plot_interpolated_n_analytical_soln_along_line(nodeCoords,elem,...
                                        edge_node,UUR,line,npts,scale,specs,isAnalytical,soln)
nAnalytical=9; % the position of the analytical function in the arguments 
if(isempty(scale))
    scale.length=1;
    scale.uo=0;
    scale.ud=1;
end
nargSpecs = 8; % position of specs
if(nargin<nargSpecs)
    specs.line = 'b-';
    specs.marker = 'bo';
    specs.width = 2;
    specs.mksize = 8;
    specs.analytical = 'r-.';
end
ftsize=20;

[xI,yI,sI,bm] = subdivide_line_by_element_edges(line,npts,nodeCoords,edge_node);

[u_duI,xI,yI,del]=interpolate_soln(nodeCoords,elem,UUR,xI,yI);
sI(del) = [];
bm(del) = [];

% find derivative along line
unitTan = line(3:4)'/norm(line(3:4),2);
duI = u_duI(:,2:3)*unitTan;


% scale quantities
u_duI(:,1) = (u_duI(:,1)-scale.uo)/scale.ud;
duI = duI/scale.ud*scale.length;
sI = sI/scale.length;
sIa = sI;

if(nargin>=nAnalytical)
    if(isAnalytical)        
        u_du = soln(xI,yI);
        du = u_du(:,2:3)*unitTan;
     % conforming mesh solution   
    else
        [u_du,~,~,del]=soln(xI,yI);        
        sIa(del) = [];
        du = u_du(:,2:3)*unitTan;
    end
    u_du(:,1) = (u_du(:,1)-scale.uo)/scale.ud;
    du = du/scale.ud*scale.length;
end


varargout{1} = [];
varargout{2} = [];
varargout{3} = [];
varargout{4} = [];
varargout{5} = [];
varargout{6} = [];

varargout{1} = plot(sI,u_duI(:,1),specs.line,'linewidth',specs.width);
hold on
varargout{2} = plot(sI(bm),u_duI(bm,1),specs.marker,'linewidth',specs.width,'markersize',specs.mksize);
if(nargin>=nAnalytical)        
    varargout{3} = plot(sIa,u_du(:,1),specs.analytical,'linewidth',specs.width);
end
xlabel('s/L','fontsize',ftsize)
ylabel('T/T_{max}','fontsize',ftsize)
set(gca,'fontsize',ftsize)
xlim([min(sI),max(sI)]);

if(isempty(varargout{2}) && isempty(varargout{3}))
    legend(varargout{1},'IGFEM')
elseif(isempty(varargout{3}))
    legend(varargout{2},'IGFEM')
else
    legend([varargout{2},varargout{3}],'IGFEM','analytical')
end
% derivative
%{
figure    
varargout{4} = plot(sI,duI,specs.line,'linewidth',specs.width);
hold on
varargout{5} = plot(sI(bm),duI(bm),specs.marker,'linewidth',specs.width);
if(nargin>=nAnalytical)
    varargout{6} = plot(sIa,du,specs.analytical,'linewidth',specs.width);
end

xlabel('s/L','fontsize',ftsize)
ylabel('L/T_{max}dT/ds','fontsize',ftsize)
set(gca,'fontsize',ftsize)
xlim([min(sI),max(sI)]);
   
if(isempty(varargout{5}) && isempty(varargout{6}))
    legend(varargout{4},'IGFEM')
elseif(isempty(varargout{6}))
    legend(varargout{5},'IGFEM')
else
    legend([varargout{5},varargout{6}],'IGFEM','analytical')
end
%}

end