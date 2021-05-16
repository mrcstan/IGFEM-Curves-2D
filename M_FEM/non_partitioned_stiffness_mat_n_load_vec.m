function [K,P] = non_partitioned_stiffness_mat_n_load_vec(Kff,...
                                                          Kfp,...
                                                          Kpf,...
                                                          Kpp,...
                                                          Pf,...
                                                          Pp,...
                                                          eqNum)
nNodes = numel(eqNum);
K = zeros(nNodes,nNodes);
P = zeros(nNodes,1);
Kff = full(Kff);
Kfp = full(Kfp);
Kpf = full(Kpf);
Kpp = full(Kpp);
fGloInd = (eqNum > 0);
pGloInd = ~fGloInd;
P(fGloInd) = Pf(eqNum(fGloInd));
P(pGloInd) = Pp(-eqNum(pGloInd));
K(fGloInd,fGloInd) = Kff(eqNum(fGloInd),eqNum(fGloInd));
K(fGloInd,pGloInd) = Kfp(eqNum(fGloInd),-eqNum(pGloInd));
K(pGloInd,fGloInd) = Kpf(-eqNum(pGloInd),eqNum(fGloInd));
K(pGloInd,pGloInd) = Kpp(-eqNum(pGloInd),-eqNum(pGloInd));
end