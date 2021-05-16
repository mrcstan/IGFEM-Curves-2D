clear all
close all
path(path, '../channelFiles')
path(path, '../M_channels')
path(path, '../../NURBS/nurbs_toolbox')
path(path, '../M_geom_toolbox')
channelFile = 'case8_parallel_stephen.channel';
%readOptions.boundsFile = 'parallel2_start_bounds.channel';
readOptions.boundsFile = [];
channels = preprocess_channels(channelFile,readOptions);

tol = 1e-10;
flag = channels_self_intersections(channels,tol)