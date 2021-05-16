close all
clear all
path(path,'../channelFiles')
path(path, '../M_geom_toolbox')
path(path, '../../NURBS/nurbs_toolbox')
path(path,'../../Opt-IGFEM-Curves-2D/M_optimization')

channelFile = 'parallel7_start_bounds_w_diams.channel';
polygonFile = 'parallelSeven_w_diams_NE.polygon';
nSamples = 24;
maxTrials = 1000;
nlconfun = @nonlinear_constraints;
%nlconfun = [];
nlcon.minPolyArea = 0.001*0.15*0.2;
nlcon.sinMinPolyAngle = sin(0.5*pi/180);
nlcon.minSidePolyArea = 0.001*0.15*0.2;
nlcon.sinMinSidePolyAngle = sin(0.5*pi/180);

[filename,pathname] = uiputfile('*.lhs','Save as');
outFile = [pathname,filename];
lhs = design_param_latin_hypercube_sampling(outFile,...
                                            channelFile,...
                                            polygonFile,...
                                            nSamples, ...
                                            maxTrials,...
                                            nlconfun,...
                                            nlcon);
                                        

                                        
