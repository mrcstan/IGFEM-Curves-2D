%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/1/2014
%%% Last modified date: 10/14/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this main function is for testing the function
% network_pressure_mass_flow_rate function, which output the pressures at
% the nodes (end points of channels) and the mass flow rates in the
% channels. Uncomment the examples that you want to see. 
clear all
close all
g = 9.81; % gravitational accelaration in m/s^2 


% CASE I
% This example illustrates the use of type 3 input of argument 2 of the
% function network_pressure_mass_flow_rate. In this case, the positions of
% the nodes are unknown. Only the channel lengths are known. 
% There are 4 nodes and 5 channels. This example is for real pipe system
% where the pipes are kilometers long. So do not be surprised to see the
% values of the pipe lenghts
%{
chan_nodes = [1,2;
              1,3;
              2,3;
              2,4;
              3,4]; % connectivity of the nodes. each row represents a channel.
                    % so row 1 is channel 1 with nodes 1 and 2 and row 5 is
                    % channel 5 with nodes 3 and 4
geom = [1000,1000,2000,2000,2000]'; % channel lengths in m (in column vector format)
diams = [0.4,0.2,0.283,0.283,0.573]'; % channel diameters in m (in column vector format)
heights = [0.4,0.2,0.283,0.283,0.573]';
crossSection = 'rectangular';
nu = 1e-6; % kinematic viscosity in m^2/s
nodalSources = [2,3]'; % nodes that are sources/sinks (in column vector format).
                       % in this case, nodes 2 and 3 have sources. in our
                       % battery panel simulation, the flow rate at the inlet 
                       % of the channels can be specified here
sourceStrengths = [10,10]'; % volume flow rates (in column vector format)
                            % at the source/sink nodes. In this case,
                            % nodes 2 and 3 both have volume flow rates of
                            % 10m^3/s
nodalBCs = [1,4]'; % nodes that have their pressures fixed (in column vector format)
density = 1000; % in kg/m^3
fixedPressures = [20,10]'*g*density; % pressures (in column vector format) in Pa

% the function network_pressure_mass_flow_rate returns
% the pressure at each node and the flow rates in each channel in column
% vector format
[pressures,flowrates] = network_pressure_mass_flow_rate(chan_nodes,...
                                                  geom,...
                                                  diams,...
                                                  heights,...
                                                  nu,...
                                                  nodalSources,...
                                                  sourceStrengths,...
                                                  powerXdensity,...
                                                  nodalBCs,...
                                                  BCs,...
                                                  crossSection);
                                              
                                                                                             

geom = [1,1;
        1,0;
        0,1;
        0,0]; % to plot the channel network, we provide fictitious coordinates
              % for the nodal locations. each row corresponds to a node.
              % for example, node 1 has coordinates (1,1), node 2 has
              % coordinates (1,0)
% the plotted network of channels will show the node number in angular brackets 
% and the pressure at the right side of the node number in blue.
% the channel number is parentheses and it is followed by the mass flow
% rate in red. 
% the arrows show the direction of the flow
plot_channel_network(chan_nodes,geom,pressures,flowrates)
%}
% CASE 2
% this example uses NURBS for the second input of network_pressure_mass_flow_rate.
% So you may skip this.
%{
p = [0,0;
     0.5,2;
     1,0]'*430.344;
knot = [0,0,0,1,1,1];
nurbs(1) = nrbmak(p,knot);
p = [0,0;
     0.5,-2;
     1,0]'*430.344;
nurbs(2) = nrbmak(p,knot);
chan_nodes = [1,2;
              1,2];
diam = [0.1,0.2]';
nu = 1e-6;
nodalSources = 2;
sourceStrengths = -2;
nodalBCs = 1;
fixedghead = 10*g;
[ghead,flowrate,len] = network_pressure_mass_flow_rate(chan_nodes,...
                                                   nurbs,...
                                                   diams,...
                                                   heights
                                                   nu,...
                                                   nodalSources,...
                                                   sourceStrengths,...
                                                   [],...
                                                   nodalBCs,...
                                                   fixedghead,...
                                                   crossSection);
head = ghead/g;
plot_channel_network(chan_nodes,nurbs,head,flowrate)
plot_channel_network(chan_nodes,nurbs,[],len)
%}
% CASE 3
% this example illustrates the use of the positions of the nodes for the
% 2nd input of network_pressure_mass_flow_rate
%{
chan_nodes = [1,2;
              2,4;
              2,5;
              2,3;
              4,5;
              3,5];          
geom = [0,0;
        500,500;
        1100,500;
        500,1200;
        1100,1200]; % location of the nodes. each row represents a node
                    % and the corresponding (x,y) coordinates
                    % for example, node 1 has coordinates (0,0)
diam = [0.15,0.1,0.075,0.1,0.15,0.15]';
nu = 9.6236e-6;
nodalSources = [3,4,5]';
sourceStrengths = [-0.2,-0.5,-0.25]'; % volume flow rates
nodalBC = 1;
density = 1000;
fixedPressure = 100*g*density; 
% if the 3rd output argument is requested, the function returns the length
% of each channel as well. In this case, this is assigned to len.
[pressures,flowrates,len] = network_pressure_mass_flow_rate(chan_nodes,...
                                                  geom,...
                                                  diams,...
                                                  heights,...
                                                  nu,...
                                                  nodalSources,...
                                                  sourceStrengths,...
                                                  [],...
                                                  nodalBC,...
                                                  fixedPressure,...
                                                  crossSection);
                                                                                            
% this time the real positions of the nodes are plotted
plot_channel_network(chan_nodes,geom,pressures,flowrates)
%}
% CASE 4
% network of parallel channels for a battery cooling panel
%{
chan_nodes = [1,2;    
              2,3;
              4,5;
              6,7;
              8,9;
              10,11;
              12,13;
              14,15;
              16,17;
              18,19;
              2,4;
              4,6;
              6,8;
              8,10;
              10,12;
              12,14;
              14,16;
              16,18;
              3,5;
              5,7;
              7,9;
              9,11;
              11,13;
              13,15;
              15,17;
              17,19;
              19,20];
               
geom =   [0,    0.189;
          0.01, 0.189;
          0.139, 0.189;
          0.01, 0.16675;
          0.139,0.16675;
          0.01, 0.1445;
          0.139,0.1445;
          0.01, 0.12225;
          0.139,0.12225;
          0.01, 0.1;
          0.139,0.1;
          0.01, 0.07775;
          0.139,0.07775;
          0.01, 0.0555;
          0.139,0.0555;
          0.01, 0.03325;
          0.139,0.03325;
          0.01, 0.011;
          0.139,0.011;
          0.15, 0.011];

diam = 7.5e-4; % if diameters are all the same, you can just specify a scalar
      
nu = 3.405e-6; % kinematic viscosity

inletNode = 1; % inlet node

inMass = 5e-4; % inlet mass flow rate kg/s

outletNode = 20; % outlet node

outPressure = 0; % outlet pressure Pass
[pressures,flowrates,len] = network_pressure_mass_flow_rate(chan_nodes,...
                                                            geom,...
                                                            diams,...
                                                            heights,...
                                                            nu,...
                                                            inletNode,...
                                                            inMass,...
                                                            [],...
                                                            outletNode,...
                                                            outPressure,...
                                                            crossSection);
% the flow rates are multiplied by 1000 to show them in g/s     
plot_channel_network(chan_nodes,geom,pressures,flowrates*1000) 
%}

% CASE 4
path(path,'../ChannelFiles')
path(path,'../M_channels')
path(path, '../../NURBS/nurbs_toolbox')
%hannelFile = 'parallel2x2_test_sensitivity.channel';
channelFile = 'parallel4_start.channel';
channels = read_channels(channelFile);


%channels.powerXdensity = 1.535225948603977e+01;
%channels.powerXdensity  = 5;
[pressure,mass] ...
        = network_pressure_mass_flow_rate(channels.contvty,...
                                          channels.nurbs,...
                                          channels.diams,...
                                          channels.heights,...
                                          channels.viscosity,...
                                          channels.inletEndPoint,...
                                          channels.massin,...
                                          channels.powerXdensity,...
                                          channels.pressureOutletEndPoint,...
                                          channels.pressureOut,...
                                          channels.crossSection);


plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)

