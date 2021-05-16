clear all
close all
path(path,'../M_Channels')

vertices = [0,0.19;
            0.15,0.01;
            0.01,0.19;
            0.14,0.19;
            0.01,0.1;
            0.14,0.1;
            0.01,0.01;
            0.14,0.01];
adjMat = generate_adj_mat_from_degree_sequence(degSeq);
gplot(adjMat,vertices,'k-','linewidth',2);