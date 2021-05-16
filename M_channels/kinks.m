%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/11/2014
%%% Last modified date: 12/25/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kinks = kinks(nurbs)
kinks = struct('x',[],'y',[],'nurbsParam',[],'seg',[],'ctrlPtNum',[]);
for i = 1:length(nurbs)
    if (nurbs(i).order == 2 && nurbs(i).number > 2)
        kinks.x = [kinks.x, nurbs(i).coefs(1,2:end-1)./nurbs(i).coefs(4,2:end-1)];
        kinks.y = [kinks.y, nurbs(i).coefs(2,2:end-1)./nurbs(i).coefs(4,2:end-1)];
        kinks.nurbsParam = [kinks.nurbsParam, ...
                            nurbs(i).knots((nurbs(i).order+1):nurbs(i).number)];
        kinks.seg = [kinks.seg, i*ones(1,nurbs(i).number-2)]; 
        kinks.ctrlPtNum = [kinks.ctrlPtNum, 2:(nurbs(i).number-1)];
    end
end
end