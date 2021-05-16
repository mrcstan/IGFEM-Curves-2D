%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/26/2013
%%% Last modified date: 9/9/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: c is the c-th child element.
% OUTPUT:
%   element-wise quantity
%   varargout(1): (L2 error)^2
%            (2): (L2 norm of analytical solution)^2 
%            (3): (H2 error)^2
%            (4): (H2 norm of analytical solution)^2
function varargout=error_squared_L2_H1_child_element(analytical,UUR,...
                                             nodeCoord,parent,gauss,c)     
varargout{1}=0;
varargout{2}=0;
if(nargout>2)
    varargout{3}=0;
    varargout{4}=0;
end


% indices of original nodes in parent element
iPON=1:3;
% indices of enrichment nodes in parent element
iPEN=4:parent.nNode;
xOriginal=[nodeCoord(parent.node(iPON),1),...
           nodeCoord(parent.node(iPON),2)]';
       
% iCEN refers to the indices of the enrichment nodes of the child
% element 
[~,iCENP,iCEN]=intersect(parent.node(iPEN),parent.child(c).node,'stable');
iCENP=iPEN(iCENP);
% iPEN are the indices of the enrichment nodes in the parent element
nurbsSurf=parent.child(c).nurbsSurf;
nuknot=length(nurbsSurf.knots{1});
nvknot=length(nurbsSurf.knots{2});
dnurbsSurf=parent.child(c).dnurbsSurf;
nurbsInd=parent.child(c).nurbsInd(iCEN);
% make nurbs shape functions
 if(~parent.child(c).isTriangle && length(parent.child(c).nurbsSeg)>1 ...
        && length(unique(parent.child(c).lineSource))==3)      
    [~,iSCEN] = ismember(parent.child(c).node(iCEN),...
                        [parent.child(c).lineSource(1),...
                        parent.child(c).lineSource(end)]);
   [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd,iSCEN);                 
else    
    [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd);
end
for i=nurbsSurf.order(1):nuknot-nurbsSurf.order(1)
    for j=nurbsSurf.order(2):nvknot-nurbsSurf.order(2)
         if(nurbsSurf.knots{1}(i)~=nurbsSurf.knots{1}(i+1) ...
           && nurbsSurf.knots{2}(j)~=nurbsSurf.knots{2}(j+1))
            rSubEl=[nurbsSurf.knots{1}(i),nurbsSurf.knots{1}(i+1),...
                    nurbsSurf.knots{1}(i+1),nurbsSurf.knots{1}(i);...
                    nurbsSurf.knots{2}(j),nurbsSurf.knots{2}(j),...
                    nurbsSurf.knots{2}(j+1),nurbsSurf.knots{2}(j+1)];
            for k = 1:gauss.qua1Dt.npt
                for m = 1:gauss.qua1Dn.npt 
                    r = [gauss.qua1Dt.pt(k),gauss.qua1Dn.pt(m)];
                    w = gauss.qua1Dt.weight(k)*gauss.qua1Dn.weight(m);
                    % calculate quantities from parametric space to natural
                    % coordinate space
                    [Np2n,DNp2n]=shape_funct(r,2);
                    Jp2n=rSubEl*DNp2n;
                    detJp2n=det(Jp2n);
                    rParam=rSubEl*Np2n;
                    % calculate quantities from physical space to parametric space              
                    % compute shape functions
                    [detJp2p,Np2p,Bp2p,~]=compute_child_element_JNB(...
                                   nurbsSurf,dnurbsSurf,nurbsShape,dnurbsShape,...
                                   iPON,iCENP,parent.nNode,xOriginal,rParam);   
                    %if(detJp2p<=0)
                    %        warning('jacobian non-positive')
                    %end              
                    factor=detJp2n*detJp2p*w;
                    XX=nrbeval(nurbsSurf,[rParam(1);rParam(2)]);
                    u_du=analytical(XX(1),XX(2));
                     if(isempty(u_du))
                         warning('integration point cannot be found in any element')
                    end
                    UU=UUR(parent.node)'*Np2p;
                    dU=UUR(parent.node)'*Bp2p;
                    UU_ua_sq=(UU-u_du(1))^2;
                    varargout{1}=varargout{1}+UU_ua_sq*factor;
                    ua_sq=u_du(1)^2;
                    varargout{2}=varargout{2}+ua_sq*factor;
                    if(nargout>2)
                        vector=[dU(1)-u_du(2),dU(2)-u_du(3)];                
                        varargout{3}=varargout{3}+(UU_ua_sq+dot(vector,vector))*factor;
                        vector=[u_du(2),u_du(3)];
                        varargout{4}=varargout{4}+(ua_sq+dot(vector,vector))*factor;
                    end
                end
            end      
         end
    end
end

end
