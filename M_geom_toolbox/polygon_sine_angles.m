%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/30/2014
%%% Last modified date: 11/30/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the angles at the polygon vertices assuming that
% the vertices are ordered CCW as well as the derivatives of the angles wrt
% the coordinates of the points
function [sineAngles,dSineAngles] = polygon_sine_angles(X,Y)
nv = numel(X);
if (nv ~= numel(Y))
    error('X and Y must have same number of elements')
end

%{
if (nargin < 3 )
    tolsq = (1e-13)^2;
else
    tolsq = tol^2;
end
%}
sineAngles = zeros(nv,1);
if (nargout > 1)
    dSineAngles = zeros(2*nv,nv);
end
ind = [nv,1:nv,1]; % wrapped index
for i = 2:nv+1
    u = [X(ind(i+1))-X(ind(i)),Y(ind(i+1))-Y(ind(i))];
    v = [X(ind(i-1))-X(ind(i)),Y(ind(i-1))-Y(ind(i))];
    normusq = dot(u,u);
    normvsq = dot(v,v);
    %if (normusq > tolsq && normvsq > tolsq)
        ucrossv = cross([u,0],[v,0]);
        normu = sqrt(normusq);
        normv = sqrt(normvsq);
        sineAngles(ind(i)) = ucrossv(3)/(normu*normv);
        if (nargout > 1)
            numer = (X(ind(i+1))-X(ind(i)))*(Y(ind(i-1))-Y(ind(i)))...
                     -(Y(ind(i+1))-Y(ind(i)))*(X(ind(i-1))-X(ind(i)));
            inv_normunormv = 1./(normu*normv);
            inv_normu3normv = 1./(normu^3*normv);
            inv_normunormv3 = 1./(normu*normv^3);

            dSineAngles(2*ind(i)-1,ind(i)) = inv_normunormv*(Y(ind(i+1))-Y(ind(i-1)))...
                                            +numer*(inv_normu3normv*u(1) ...
                                                   +inv_normunormv3*v(1));

            dSineAngles(2*ind(i),ind(i)) = inv_normunormv*(X(ind(i-1))-X(ind(i+1)))...
                                          +numer*(inv_normu3normv*u(2) ...
                                                 +inv_normunormv3*v(2));                      
            
            dSineAngles(2*ind(i+1)-1,ind(i)) = inv_normunormv*v(2)...
                                              -numer*inv_normu3normv*u(1);
            dSineAngles(2*ind(i+1),ind(i)) = -inv_normunormv*v(1) ...
                                             -numer*inv_normu3normv*u(2);
            
            dSineAngles(2*ind(i-1)-1,ind(i)) = -inv_normunormv*u(2)...
                                               -numer*inv_normunormv3*v(1);
            dSineAngles(2*ind(i-1),ind(i)) = inv_normunormv*u(1)...
                                              -numer*inv_normunormv3*v(2);
        end
    %end
end

%{
% i = 1
u = [X(2)-X(1),Y(2)-Y(1)];
v = [X(nv)-X(1),Y(nv)-Y(1)];
normu = norm(u);
normv = norm(v);
if (normu > tol && normv > tol)
    ucrossv = cross([u,0],[v,0]);
    sineAngles(1) = ucrossv(3)/(normu*normv);
end
% i = nv
u = [X(1)-X(nv),Y(1)-Y(nv)];
v = [X(nv-1)-X(nv),Y(nv-1)-Y(nv)];
normu = norm(u);
normv = norm(v);
if (normu > tol && normv > tol)
    ucrossv = cross([u,0],[v,0]);
    sineAngles(nv) = ucrossv(3)/(normu*normv);
end
%}
end