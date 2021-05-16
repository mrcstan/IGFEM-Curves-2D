scale.length = 0.1;
scale.uo = 20.0;
scale.ud = 2000*0.1/(2*0.6);
npts = 1;


specs.mksize = 8.0;
specs.width = 2;
        
line = [0,0.0255,1,0];


%UUR2 = update_enrichment_node_value(UURA,meshA.node,meshA.edge,meshA.elem);   
%outfile = '10x5.vtk';
%scalarname = 'T';
%matlab2vtk_scalar(outfile,scalarname,meshA.elem,meshA.node,UUR2);
figure
tic
specs.line = 'b-';
specs.marker = 'bo';
[~,h] =  plot_interpolated_n_analytical_soln_along_line(meshA.node.coords,...
    meshA.elem,meshA.edge,UURA,line,npts*5,scale,specs);     
toc            
hold on

tic
specs.line = 'g-';
specs.marker = 'g^';
[~,hout] =  plot_interpolated_n_analytical_soln_along_line(meshB.node.coords,...
    meshB.elem,meshB.edge,UURB,line,npts*5,scale,specs);     
toc  
h = [h,hout];

%
tic
specs.line = 'k-';
specs.marker = 'ks';
[~,hout] =  plot_interpolated_n_analytical_soln_along_line(meshC.node.coords,...
    meshC.elem,meshC.edge,UURC,line,npts*5,scale,specs);     
toc  
h = [h,hout];
%
%
tic
specs.line = 'm-';
specs.marker = 'mv';
[~,hout] =  plot_interpolated_n_analytical_soln_along_line(meshD.node.coords,...
    meshD.elem,meshD.edge,UURD,line,npts,scale,specs);     
toc  
h = [h,hout];
%
%
tic
specs.line = 'c-';
specs.marker = 'c<';
[~,hout] =  plot_interpolated_n_analytical_soln_along_line(meshE.node.coords,...
    meshE.elem,meshE.edge,UURE,line,npts*5,scale,specs);     
toc  
h = [h,hout];
%
%
tic
specs.line = 'y-';
specs.marker = 'y>';
[~,hout] =  plot_interpolated_n_analytical_soln_along_line(meshF.node.coords,...
    meshF.elem,meshF.edge,UURF,line,npts,scale,specs);     
toc  
h = [h,hout];
%
%{
tic
specs.line = 'r-';
specs.marker = 'rd';
[hout,~] =  plot_interpolated_n_analytical_soln_along_line(meshFine.node.coords,...
    meshFine.elem,meshFine.edge,UURFine,line,npts,scale,specs); 
toc                            
h = [h,hout];
%}                     
%legend(h,'IGFEM (n_{el}=100)','SFEM (n_{el}=116)','IGFEM (n_{el}=400)','SFEM (n_{el}=416)','SFEM (n_{el}=26944)')                         
%
legend(h,'IGFEM (n_{el}=100)','SFEM (n_{el}=122)','IGFEM (n_{el}=400)','SFEM (n_{el}=458)','IGFEM (n_{el}=1600)','SFEM (n_{el}=1832)','SFEM (n_{el}=26944)')
%

