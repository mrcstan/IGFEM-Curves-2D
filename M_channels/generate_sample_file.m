close all
clear all
path(path,'../ChannelFiles')
path(path,'../GeomConstrFiles')
path(path, '../M_geom_toolbox')
path(path, '../../NURBS/nurbs_toolbox')
path(path,'../../Opt-IGFEM-2D/M_optimization')
path(path,'../../vflib/mex')
path(path,'../../MatlabUsefulFunctions/GrTheory')

samplingType = 'random0';
% must provide non-overlapping bounding boxes for samplingType = 'random0' and 'lhs' 
channelFile = '3x3ref.channel'; 
polygonFile = '3x3NE.polygon';
nSamples = 24;
maxTrials = 3000;
%nlconfun = @nonlinear_constraints;
nlconfun = [];
nlcon.minPolyArea = 0.001*0.15*0.2;
nlcon.sinMinPolyAngle = sin(0.5*pi/180);
nlcon.minSidePolyArea = 0.001*0.15*0.2;
nlcon.sinMinSidePolyAngle = sin(0.5*pi/180);
rng('shuffle')    

% needed for both 'random1' and 'subgraph_monomorphism'
coordBounds = [0.005,0.145;0.005,0.195]; 

% only valid for samplingType = 'random1';
nBranches = 8; % REMEMBER to choose correct number of branches for samplingType = 'random1'
length1 = 0.08/nBranches;
length2 = 0.12/nBranches;
lengthBounds = [0.1,0.14;
                length1,length2;
                0.1,0.14;
                length1,length2];
minAngle = 5*pi/180;
if (strcmpi(samplingType,'random0'))
    [filename,pathname] = uiputfile('*.rand','Save as','../SampleFiles/');
    outFile = [pathname,filename];
    sample = design_param_random_sampling(channelFile,...
                                          polygonFile,...
                                          nSamples, ...
                                          maxTrials,...
                                          nlconfun,...
                                          nlcon);
elseif (strcmpi(samplingType,'random1'))
    [filename,pathname] = uiputfile('*.rand1','Save as','../SampleFiles/');
    outFile = [pathname,filename];
  
    sample = design_param_random_parallel_network(channelFile, ...
                                                  polygonFile,...
                                                  nSamples, ...
                                                  maxTrials, ...
                                                  coordBounds, ...
                                                  lengthBounds, ...
                                                  minAngle,...
                                                  nlconfun,...
                                                  nlcon);
                                            
elseif (strcmpi(samplingType,'random2'))   
    [filename,pathname] = uiputfile('*.rand2','Save as','../SampleFiles/');
    outFile = [pathname,filename];
    sample = design_param_random_assign_ctrl_pts(channelFile,...
                                                  polygonFile,...
                                                  nSamples, ...
                                                  maxTrials,...
                                                  nlconfun,...
                                                  nlcon);
elseif (strcmpi(samplingType,'subgraph_monomorphism'))
    [filename,pathname] = uiputfile('*.smmp','Save as','../SampleFiles/');
    outFile = [pathname,filename];
    %outFile = 'test.smmp';
    sample = design_param_subgraph_monomorphism(channelFile, ...
                                                  polygonFile,...
                                                  nSamples, ...
                                                  maxTrials, ...
                                                  coordBounds, ...
                                                  nlconfun,...
                                                  nlcon);                                              
elseif (strcmpi(samplingType,'lhs'))
    [filename,pathname] = uiputfile('*.lhs','Save as','../SampleFiles/');
    outFile = [pathname,filename];
    sample = design_param_latin_hypercube_sampling(channelFile,...
                                                   polygonFile,...
                                                   nSamples, ...
                                                   maxTrials,...
                                                   nlconfun,...
                                                   nlcon);
else
    error('unknown sampling type')
end
                                 
nSamples = size(sample,2);
nParams = size(sample,1);
fileID = fopen(outFile,'w');
fprintf(fileID,'nSamples,%i\n',nSamples);
fprintf(fileID,'nParams,%i\n',nParams);
for j = 1:nSamples
    fprintf(fileID,'%g \n',sample(:,j));
end
fclose(fileID);

