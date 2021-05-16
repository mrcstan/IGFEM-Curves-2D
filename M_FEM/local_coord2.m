%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan on 7/17/2013
%%% Last modified date: 07/17/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function returns local coordinates
% INPUT:
% shape: 1 - triangular element
%        2 - quadrilateral element
% NOTE:
% two methods to obtain local coordinates of quadrilateral element
% i) directly solving the equations (not fully debugged, use with caution)
% ii) using fsolve (preferred). only involves the function
% invers_mapping_quadrilateral.

function Xloc = local_coord(shape, Xglo, Xel)
Xloc=zeros(2,1);
if shape == 1
    A=[Xel(1,2)-Xel(1,1),Xel(1,3)-Xel(1,1);Xel(2,2)-Xel(2,1),Xel(2,3)-Xel(2,1)]; 
    rhs = [Xglo(1)-Xel(1,1); Xglo(2)-Xel(2,1)];
    Xloc=A\rhs;    
elseif shape == 2
    %{
    %Xglo(1)=a0+a1*Xloc(1)+a2*Xloc(2)+a3*Xloc(1)*Xloc(2)
    %Xglo(2)=b0+b1*Xloc(1)+b2*Xloc(2)+b3*Xloc(1)*Xloc(2)
    a0=sum(Xel(1,1:4));
    a1=-Xel(1,1)+Xel(1,2)+Xel(1,3)-Xel(1,4);
    a2=-Xel(1,1)-Xel(1,2)+Xel(1,3)+Xel(1,4);
    a3=Xel(1,1)-Xel(1,2)+Xel(1,3)-Xel(1,4);
    b0=sum(Xel(2,1:4));
    b1=-Xel(2,1)+Xel(2,2)+Xel(2,3)-Xel(2,4);
    b2=-Xel(2,1)-Xel(2,2)+Xel(2,3)+Xel(2,4);
    b3=Xel(2,1)-Xel(2,2)+Xel(2,3)-Xel(2,4);
    if(a3==0 && b3==0)
        A=[a1,a2;b1,b2];
        rhs=[4*Xglo(1)-a0;4*Xglo(2)-b0];
        Xloc=A\rhs;
    else
        c0=a0*b3-a3*b0-4*(Xglo(1)*b3-Xglo(2)*a3);
        c1=a1*b3-a3*b1;
        c2=a2*b3-a3*b2;
        if(abs(c2)>tol)
            d2=-a3*c1/c2;
            d1=a1-a2*c1/c2-a3*c0/c2;
            d0=a0-a2*c0/c2-4*Xglo(1);
            if(abs(d2)<tol)
                d2=0.0;
            end
            if(abs(d0)<tol && abs(d1)<tol && abs(d2)<tol)
                d2=-b3*c1/c2;
                d1=b1-b2*c1/c2-b3*c0/c2;
                d0=b0-b2*c0/c2-4*Xglo(2);
                flag=false;
            end
            rt=roots([d2,d1,d0]);
            for i=1:length(rt)
               xi=rt(i);
               eta=-(c0+c1*xi)/c2;
               if(flag)
                   res=4*Xglo(2)-b0-b1*xi-b2*eta-b3*xi*eta;
                   if(abs(res)<tol2);
                        Xloc(1)=xi;
                        Xloc(2)=eta;
                        return
                   end
               else
                   res=4*Xglo(1)-a0-a1*xi-a2*eta-a3*xi*eta;
                   if(abs(res)<tol2);
                        Xloc(1)=xi;
                        Xloc(2)=eta;
                        return
                   end
               end
            end
        elseif(abs(c1)>tol)
            d2=-a3*c2/c1;
            d1=a2-a1*c2/c1-a3*c0/c1;
            d0=a0-a1*c0/c1-4*Xglo(1);
            if(abs(d0)<tol && abs(d1)<tol && abs(d2)<tol)
                d2=-b3*c2/c1;
                d1=b2-b1*c2/c1-b3*c0/c1;
                d0=b0-b1*c0/c1-4*Xglo(2);
                flag=false;
            end
            if(abs(d2)<tol)
                d2=0.0;
            end
            rt=roots([d2,d1,d0]);
            for i=1:length(rt)
               eta=rt(i);
               xi=-(c0+c2*eta)/c1;
               if(flag)
                   res=4*Xglo(2)-b0-b1*xi-b2*eta-b3*xi*eta;
                   if(abs(res)<tol);
                        Xloc(1)=xi;
                        Xloc(2)=eta;
                        return
                   end
               else
                   res=4*Xglo(1)-a0-a1*xi-a2*eta-a3*xi*eta;
                   if(abs(res)<tol);
                        Xloc(1)=xi;
                        Xloc(2)=eta;
                        return
                   end
               end
            end
        else
            error('local_coord: either c1 or c2 must not vanish')
        end
    end
    %}
    Xloc = invers_mapping_quadrilateral(Xglo, Xel);
end

end

