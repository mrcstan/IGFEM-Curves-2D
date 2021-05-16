scale.length = 0.1;
scale.uo = 0.0;
scale.ud = 1.0;
npts = 501;

xi = min(mesh.node.coords(:,1));
xf = max(mesh.node.coords(:,1));
yi = min(mesh.node.coords(:,2));
yf = max(mesh.node.coords(:,2));
%
xI=linspace(xi,xf,npts)';
yI=0.01*ones(npts,1);

figure
xMark = linspace(xi,xf,11);
specs.line = 'b-';
specs.marker = 'bo';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshA.node.coords,meshA.elem,UURA,...
                                xI,yI,'x',xMark,scale,specs);
toc
hold on
h = hout;

%
xMark = linspace(xi,xf,11);
specs.line = 'g-';
specs.marker = 'g^';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshB.node.coords,meshB.elem,UURB,...
                                xI,yI,'x',xMark,scale,specs);       
toc                                                     
h = [h,hout];

%
xMark = linspace(xi,xf,21);
specs.line = 'k-';
specs.marker = 'ks';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshC.node.coords,meshC.elem,UURC,...
                                xI,yI,'x',xMark,scale,specs);    
toc                                                    
h = [h,hout];

xMark = linspace(xi,xf,21);
specs.line = 'm-';
specs.marker = 'mv';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshD.node.coords,meshD.elem,UURD,...
                                xI,yI,'x',xMark,scale,specs);       
toc                                                     
h = [h,hout];

xMark = linspace(xi,xf,41);
specs.line = 'c-';
specs.marker = 'c<';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshE.node.coords,meshE.elem,UURE,...
                                xI,yI,'x',xMark,scale,specs);    
toc                                                     

h = [h,hout];
xMark = linspace(xi,xf,41);
specs.line = 'y-';
specs.marker = 'y>';
tic
[~,hout] = plot_interpolated_n_analytical_soln(meshF.node.coords,meshF.elem,UURF,...
                                xI,yI,'x',xMark,scale,specs);       
toc                                                     
h = [h,hout];

%{
specs.line = 'y-';
specs.marker = 'y>';
hout = plot_interpolated_n_analytical_soln(meshFine.node.coords,meshFine.elem,UURFine,...
                                xI,yI,'x',[],scale);                            

h = [h,hout];
%}
%                            
                     
  
%
legend(h,'IGFEM (n_{el}=100)','SFEM (n_{el}=122)','IGFEM (n_{el}=400)','SFEM (n_{el}=458)','IGFEM (n_{el}=1600)','SFEM (n_{el}=1854)','SFEM (n_{el}=26676)')