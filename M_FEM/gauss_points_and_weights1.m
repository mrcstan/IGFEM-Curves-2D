%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 1/2/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns guss points and weights for triangular and
% quadrilateral elements
% isTriangular: true - the integration domain is triangular or tetrahedral 
% npt: total number of gauss points for triangular or tetrahedral
%      integration. number of gauss points in one direction for quadrilateral 
%      or hexahedral integration
function gauss=gauss_points_and_weights(isTriangular,npt,dimension)                                
    if (isTriangular)
        gauss.npt=npt;
        switch dimension
            case 1
                [gauss.pt,gauss.weight]=lgwt(npt,0,1);               
            case 2
                switch npt
                    case 3 % error is of order O(h^3)
                        gauss.pt = [1./6.,1./6.;
                                       2./3.,1./6.;
                                       1./6.,2./3.];
                        gauss.weight = [1./6., 1./6., 1./6.];
                    case 7 % error is of order O(h^6)
                        r1=0.1012865073235;
                        r2=0.7974269853531;
                        r4=0.4701420641051;
                        r6=0.0597158717898;
                        r7=1./3.;
                        gauss.pt = [r1,r1;
                                    r2,r1;
                                    r1,r2;
                                    r4,r6;
                                    r4,r4;
                                    r6,r4;
                                     r7,r7];
                        w1=0.1259391805448;
                        w4=0.1323941527885;
                        w7=0.225;
                        gauss.weight = 0.5*[w1,w1,w1,w4,w4,w4,w7];
                    otherwise
                        error('gauss_points_and_weights: number of gauss points unavailable')
                end
            case 3
                switch npt
                    case 4 % error is of order O(h^3)
                        r1 = 0.585410196624969;
                        r2 = 0.138196601125011;
                        w1 = 1./24.;
                        gauss.pt = [r1,r2,r2;
                                    r2,r1,r2;
                                    r2,r2,r1;
                                    r2,r2,r2];
                        gauss.weight = [w1,w1,w1,w1];
                    case 5 % error is of order O(h^4)
                        r1 = 0.5;
                        r2 = 1./6.;
                        w1 = -0.8/6.;
                        w2 = 0.45/6.;
                        gauss.pt = [0.25,0.25,0.25;
                                    r1,r2,r2;
                                    r2,r1,r2;
                                    r2,r2,r1;
                                    r2,r2,r2];
                        gauss.weight = [w1,w2,w2,w2,w2];                        
                    otherwise
                        error('gauss_points_and_weights: number of gauss points unavailable')    
                end  
            otherwise
                error('gauss_points_and_weights: dimension unavailable')
        end
    else     
        [pt,weight]=lgwt(npt,-1,1);
        switch dimension
            case 1
                gauss.pt=pt;
                gauss.weight=weight;
                gauss.npt=npt;
            case 2
                gauss.npt=npt*npt;
                gauss.pt(1:gauss.npt,1:2)=0;
                for i=1:npt
                    for j=1:npt
                        k=npt*(i-1)+j;
                        gauss.pt(k,1:2)=[pt(i),pt(j)];
                        gauss.weight(k)=weight(i)*weight(j);
                    end
                end                 
            case 3
                nptsq = npt*npt;
                gauss.npt=nptsq*npt;
                gauss.pt(1:gauss.npt,1:3)=0;
                for i=1:npt
                    for j=1:npt
                        for k=1:npt
                            p = nptsq*(i-1)+npt*(j-1)+k;
                            gauss.pt(p,1:3) = [pt(i),pt(j),pt(k)];
                            gauss.weight(p) = weight(i)*weight(j)*weight(k);
                        end
                    end
                end
            otherwise
                error('gauss_points_and_weights: number of gauss points unavailable')
         end
    end    
end

