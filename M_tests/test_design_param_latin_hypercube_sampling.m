close all
clear all
path(path,'../M_channels')
path(path,'../channelFiles')
channelFile = 'parallel2_start_bounds_w_diams.channel';
nSamples = 4;
maxTrials = 10;
nlconfun = @nonlinear_constraints;
nlcon.minPolyArea = 0.001*0.15*0.2;
%nlcon.minPolyArea = 1e-4;
nlcon.sinMinPolyAngle = sin(0.5*pi/180);
outFile = 'parallelTwo.lhs';
lhs = design_param_latin_hypercube_sampling(outFile,...
                                            channelFile,...
                                            nSamples, ...
                                            maxTrials,...
                                            nlconfun,...
                                            nlcon);
