path(path, './M_FEM')

% case 1
%{
Xel = [0,0;
       1,0;
       0,1]';
[~, DN] = shape_funct([0,0], 1);
J = Xel*DN; % Xel is a 2x3 matrix, DN is a 3x2 matrix
B = DN/J;
polarAngle = 90*pi/180.0;
unitTan = [cos(polarAngle);sin(polarAngle)];
he = streamwise_elem_length(unitTan,B);

hsoln = abs(sec(polarAngle)/(1+tan(polarAngle)));

relErr = abs(1-he/hsoln)
%}
% case 2
Xel = [0,0;
       0.5*sqrt(3),-0.5;
       0.5*sqrt(3),0.5]';
[~, DN] = shape_funct([0,0], 1);
J = Xel*DN; % Xel is a 2x3 matrix, DN is a 3x2 matrix
B = DN/J;
%polarAngle = 0*pi/180.0;
%unitTan = [cos(polarAngle);sin(polarAngle)];
unitTan = [0.5*sqrt(3);0.25];
unitTan = unitTan/norm(unitTan);
he = streamwise_elem_length(unitTan,B)
hsoln = sqrt(13)/4.0